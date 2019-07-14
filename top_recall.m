function top_recall
seed=1;
cv_num=5;
step=50;
CrossValidation(seed, cv_num, step)
end

function final_result=CrossValidation(seed, cv_num, step)
load('lncRNA_miRNA_all.mat');
interaction_matrix=interactionmatrix;         %lncRNA端
[row_index,col_index]=find(interaction_matrix==1);
link_num=sum(sum(interaction_matrix));
rand('state',seed);
random_index=randperm(link_num);
size_of_cv=round(link_num/cv_num);
result=zeros(1, step);
for k=1:cv_num
    if (k~=cv_num)
        test_row_index=row_index(random_index((size_of_cv*(k-1)+1):(size_of_cv*k)));
        test_col_index=col_index(random_index((size_of_cv*(k-1)+1):(size_of_cv*k)));
    else
        test_row_index=row_index(random_index((size_of_cv*(k-1)+1):end));
        test_col_index=col_index(random_index((size_of_cv*(k-1)+1):end));
    end
    train_interaction_matrix=interaction_matrix;
    test_link_num=size(test_row_index,1);
    for i=1:test_link_num
        train_interaction_matrix(test_row_index(i),test_col_index(i))=0;
    end
% %     seq
%     similairty_matrix_1=GetImprovedLNSimilarity(train_interaction_matrix, L_Kmer_5, 0.8);
%     score_matrix_1=LabelPropagation(similairty_matrix_1,train_interaction_matrix,0.4);
%     similairty_matrix_2=GetImprovedLNSimilarity(train_interaction_matrix', M_Kmer_5, 0.8);
%     score_matrix_2=LabelPropagation(similairty_matrix_2,train_interaction_matrix',0.4);
%     score_matrix = 0.25 * score_matrix_1 + (1 - 0.25) * score_matrix_2';
    
% %     exp
%     similairty_matrix_1=GetImprovedLNSimilarity(train_interaction_matrix, lncRNA_expression, 0.1);
%     score_matrix_1=LabelPropagation(similairty_matrix_1,train_interaction_matrix,0.1);
%     similairty_matrix_2=GetImprovedLNSimilarity(train_interaction_matrix', miRNA_expression, 0.1);
%     score_matrix_2=LabelPropagation(similairty_matrix_2,train_interaction_matrix',0.1);
%     score_matrix = 0.45 * score_matrix_1 + (1 - 0.45) * score_matrix_2';
    
% %     RA
%     score_matrix = ResourceAllocation(train_interaction_matrix);
    
%     CF
    similairty_matrix_1=GetImprovedCosineimilarity(train_interaction_matrix);
    score_matrix_1=CollaborativeFiltering(similairty_matrix_1,train_interaction_matrix,1);
    similairty_matrix_2=GetImprovedCosineimilarity(train_interaction_matrix');
    score_matrix_2=CollaborativeFiltering(similairty_matrix_2,train_interaction_matrix',1);
    score_matrix = 0.5 * score_matrix_1 + (1 - 0.5) * score_matrix_2';
    
    result(1,:) = result(1,:) + ModelEvaluate(interaction_matrix, score_matrix, train_interaction_matrix, step);
end
final_result=result/cv_num;
end

function similarity_matrix=GetImprovedCosineimilarity(train_interaction_matrix)
index_0_row=find(all(train_interaction_matrix==0,2)==1);

%初始化相似度矩阵
similarity_matrix=ones(size(train_interaction_matrix, 1));
for i=1:size(index_0_row, 1)
    similarity_matrix(index_0_row(i, 1), :)=zeros(1, size(train_interaction_matrix, 1));
    similarity_matrix(:, index_0_row(i, 1))=zeros(size(train_interaction_matrix, 1), 1);
end

%填充反应谱相似度
temp_interaction=train_interaction_matrix;
temp_interaction(index_0_row, :)=[];
temp_similarity_matrix=GetCosineSimilarity(temp_interaction);
count = 1;
for j=1:size(similarity_matrix(:), 1)
    if(similarity_matrix(j)==1)
        similarity_matrix(j)=temp_similarity_matrix(count);
        count=count+1;
    end
end

%填充特征相似度
% feature_similarity_matrix=GetLNSimilarity(feature_matrix, round(size(feature_matrix, 1)*neighbor_alpha));
% for k=1:size(index_0_row, 1)
%     similarity_matrix(index_0_row(k, 1), :)=feature_similarity_matrix(index_0_row(k, 1), :);
% end
end

function score_matrix=ResourceAllocation(train_interaction_matrix)
[~,miRNA_num]=size(train_interaction_matrix);
resource_allocate_matrix=zeros(miRNA_num,miRNA_num);
lncRNA_degree=sum(train_interaction_matrix,2);
miRNA_degree=sum(train_interaction_matrix,1);
for i=1:miRNA_num
    for j=1:miRNA_num
        z=0;
        set1=find(train_interaction_matrix(:,i)==1);
        set2=find(train_interaction_matrix(:,j)==1);
        set = intersect(set1,set2);
        if  ~isempty(set)
            num=size(set,1);
            for p=1:num
                if lncRNA_degree(p)~=0
                    z=z+train_interaction_matrix(set(p,1),i)*train_interaction_matrix(set(p,1),j)/lncRNA_degree(set(p,1));
                end
            end
        end
        if miRNA_degree(j)~=0
            resource_allocate_matrix(i,j)=z/miRNA_degree(j);
        end
    end
end
score_matrix=(resource_allocate_matrix*train_interaction_matrix')';
end

function similarity_matrix=GetImprovedLNSimilarity(train_interaction_matrix,feature_matrix,neighbor_alpha)
index_0_row=find(all(train_interaction_matrix==0,2)==1);

%初始化相似度矩阵
similarity_matrix=ones(size(train_interaction_matrix, 1));
for i=1:size(index_0_row, 1)
    similarity_matrix(index_0_row(i, 1), :)=zeros(1, size(train_interaction_matrix, 1));
    similarity_matrix(:, index_0_row(i, 1))=zeros(size(train_interaction_matrix, 1), 1);
end

%填充反应谱相似度
temp_interaction=train_interaction_matrix;
temp_interaction(index_0_row, :)=[];
temp_similarity_matrix=GetLNSimilarity(temp_interaction, round(size(temp_interaction, 1)*neighbor_alpha));
count = 1;
for j=1:size(similarity_matrix(:), 1)
    if(similarity_matrix(j)==1)
        similarity_matrix(j)=temp_similarity_matrix(count);
        count=count+1;
    end
end

%填充特征相似度
feature_similarity_matrix=GetLNSimilarity(feature_matrix, round(size(feature_matrix, 1)*neighbor_alpha));
for k=1:size(index_0_row, 1)
    similarity_matrix(index_0_row(k, 1), :)=feature_similarity_matrix(index_0_row(k, 1), :);
end
end

function score_matrix=LabelPropagation(W,Y,alpha)
%W==similarity_matrix; Y==train_interaction_matrix
score_matrix=(1-alpha)*pinv(eye(size(W,1))-alpha*W)*Y;
end

function result=ModelEvaluate(interaction_matrix, score_matrix, train_interaction_matrix, step)
result=zeros(1,step);
real_score=interaction_matrix(:);
predict_score=score_matrix(:);
index=train_interaction_matrix(:);
test_index=find(index==0);
real_score=real_score(test_index);
predict_score=predict_score(test_index);
sort_predict_score=sort(predict_score,'descend');
% final_num=size(real_score,1);       %按比例算
final_num=5000;       %按数量算
for i=1:step
    result(1,i)=get_recall(real_score, predict_score, sort_predict_score, int32(final_num/step*i));
end
end

function recall=get_recall(real_score,predict_score,sort_predict_score,top_number)
real_label=real_score;
predict_label=(predict_score>=sort_predict_score(top_number));

tp_index=find(real_label==1 & predict_label==1);
tp=size(tp_index,1);

tn_index=find(real_label==0 & predict_label==0);
tn=size(tn_index,1);

fp_index=find(real_label==0 & predict_label==1);
fp=size(fp_index,1);

fn_index=find(real_label==1 & predict_label==0);
fn=size(fn_index,1);

recall=tp/(tp+fn);
end

function W=GetLNSimilarity(feature_matrix,neighbor_num)
%线性邻居相似度的快速计算方法(W==similarity_matrix)
iteration_max=40;
mu=3;
X=feature_matrix;

%用欧氏距离计算邻居
row_num=size(X,1);
distance_matrix=pdist2(X,X,'euclidean');
e=ones(row_num,1);
distance_matrix=distance_matrix+diag(e*inf);
[~, si]=sort(distance_matrix,2,'ascend');
nearst_neighbor_matrix=zeros(row_num,row_num);
index=si(:,1:neighbor_num);
for i=1:row_num
    nearst_neighbor_matrix(i,index(i,:))=1;
end

%迭代求相似矩阵
C=nearst_neighbor_matrix;
rand('state',1);
W=rand(row_num,row_num);
W=(C.*W);
lamda=mu*e;
P=X*X'+lamda*e';
for i=1:iteration_max
    Q=(C.*W)*P;
    W=(C.*W).*P./Q;
    W(isnan(W))=0;
end
end

function similarity_matrix = GetCosineSimilarity(feature_matrix)
numerator_matrix = feature_matrix * feature_matrix';                %分子
distance_matrix = sum(feature_matrix .* feature_matrix, 2);
denominator_matrix = sqrt(distance_matrix * distance_matrix');          %分母
similarity_matrix = numerator_matrix ./ denominator_matrix;
similarity_matrix(isnan(similarity_matrix)) = 0;
for i = 1 : size(feature_matrix, 1)
    similarity_matrix(i, i) = 0;
end
similarity_matrix = MatrixNormalize(similarity_matrix);
end

function score_matrix = CollaborativeFiltering(similarity_matrix,train_interaction_matrix,alpha)
%前alpha的保留原值，后1-alpha的置为0
non_neighbor_num=ceil(size(similarity_matrix,1)*(1-alpha));
[~, si]=sort(similarity_matrix,2,'ascend');
index=si(:,1:non_neighbor_num);
for i=1:size(similarity_matrix,1)
    similarity_matrix(i,index(i,:))=0;
end

similarity_sum_matrix = sum(similarity_matrix, 2);
similarity_sum_diagonal_matrix = pinv(diag(similarity_sum_matrix));
row_normalized_similarity_matrix = similarity_sum_diagonal_matrix * similarity_matrix;
score_matrix = row_normalized_similarity_matrix * train_interaction_matrix;
end

function similarity_matrix=MatrixNormalize(similarity_matrix)
%把矩阵双向归一化
similarity_matrix(isnan(similarity_matrix))=0;
[row,col]=size(similarity_matrix);
for i=1:row
    similarity_matrix(i,i)=0;
end
if row==col
    for round=1:10
        d=diag(sum(similarity_matrix,2));
        d1=pinv(sqrt(d));
        similarity_matrix=d1*similarity_matrix*d1;
    end
else
    for j=1:size(similarity_matrix,1)
        if sum( similarity_matrix(j,:))~=0
            similarity_matrix(j,:) = similarity_matrix(j,:)./ sum(similarity_matrix(j,:));
        else
            similarity_matrix(j,:)=zeros(1,size(similarity_matrix,2));
        end
    end
end
end