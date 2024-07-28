# train model
nmb_cluster=200
batch=128
epochs=200
lr=0.01
data_path='expr.csv'
ckpt_path='train_res'

python main.py --nmb_cluster $nmb_cluster --batch $batch \
 --epochs $epochs --lr $lr --data_path $data_path \
 --ckpt_path $ckpt_path &