# 要素数を変えたいとき
- [main.go](input/main.go)の`elementD`
- [params.h](./ELEMENT_NUM)の`ELEMENT_NUM`

を編集する
  - TODO : これらのリソースを1元化する

# 出力する波形の時間経過を変えたい時
- [output.dat](output/dat/iteratedTime.dat)の整数を変更する

# 計算の実行（可視化まで行う）
```shell
make
```

# 可視化
```shell
python3 visualize.py
```
