OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.36800185) q[0];
sx q[0];
rz(-0.79080963) q[0];
sx q[0];
rz(0.33413449) q[0];
rz(-0.45733991) q[1];
sx q[1];
rz(-0.9442803) q[1];
sx q[1];
rz(1.2184719) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63872913) q[0];
sx q[0];
rz(-1.6067061) q[0];
sx q[0];
rz(-2.3907651) q[0];
rz(-pi) q[1];
rz(-1.6115509) q[2];
sx q[2];
rz(-2.0353122) q[2];
sx q[2];
rz(1.0711311) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.76518607) q[1];
sx q[1];
rz(-1.1131439) q[1];
sx q[1];
rz(-0.54995579) q[1];
rz(-pi) q[2];
rz(0.18913194) q[3];
sx q[3];
rz(-1.7889708) q[3];
sx q[3];
rz(1.8060341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5228287) q[2];
sx q[2];
rz(-0.4814119) q[2];
sx q[2];
rz(-0.5775601) q[2];
rz(-1.9918359) q[3];
sx q[3];
rz(-1.3883608) q[3];
sx q[3];
rz(0.66453385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98786551) q[0];
sx q[0];
rz(-2.5550714) q[0];
sx q[0];
rz(-0.38744774) q[0];
rz(0.93915835) q[1];
sx q[1];
rz(-0.99717957) q[1];
sx q[1];
rz(1.4025677) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5706022) q[0];
sx q[0];
rz(-1.5413069) q[0];
sx q[0];
rz(-1.5667314) q[0];
x q[1];
rz(2.8005373) q[2];
sx q[2];
rz(-2.555072) q[2];
sx q[2];
rz(1.3618493) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.84239292) q[1];
sx q[1];
rz(-1.0789011) q[1];
sx q[1];
rz(-2.030034) q[1];
rz(-pi) q[2];
x q[2];
rz(0.6321294) q[3];
sx q[3];
rz(-1.4940133) q[3];
sx q[3];
rz(2.0413105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.42276057) q[2];
sx q[2];
rz(-1.3964802) q[2];
sx q[2];
rz(2.823901) q[2];
rz(-0.20673949) q[3];
sx q[3];
rz(-0.59967774) q[3];
sx q[3];
rz(-0.81682214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3690255) q[0];
sx q[0];
rz(-1.7018397) q[0];
sx q[0];
rz(1.4136219) q[0];
rz(-2.6638022) q[1];
sx q[1];
rz(-1.7910035) q[1];
sx q[1];
rz(-0.40107045) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51940489) q[0];
sx q[0];
rz(-1.6739968) q[0];
sx q[0];
rz(-1.897057) q[0];
rz(0.24984078) q[2];
sx q[2];
rz(-0.65953883) q[2];
sx q[2];
rz(1.2741054) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.229278) q[1];
sx q[1];
rz(-0.73017987) q[1];
sx q[1];
rz(0.0095403949) q[1];
rz(1.1242261) q[3];
sx q[3];
rz(-2.7105769) q[3];
sx q[3];
rz(-2.9197846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5473189) q[2];
sx q[2];
rz(-1.6197562) q[2];
sx q[2];
rz(-0.55580124) q[2];
rz(2.1650971) q[3];
sx q[3];
rz(-2.591811) q[3];
sx q[3];
rz(-0.78021375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7926086) q[0];
sx q[0];
rz(-1.6225092) q[0];
sx q[0];
rz(1.4439616) q[0];
rz(-1.5199039) q[1];
sx q[1];
rz(-0.65602055) q[1];
sx q[1];
rz(-0.25340432) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4356416) q[0];
sx q[0];
rz(-2.5905529) q[0];
sx q[0];
rz(-1.7680697) q[0];
rz(-pi) q[1];
x q[1];
rz(0.01732145) q[2];
sx q[2];
rz(-2.0565196) q[2];
sx q[2];
rz(-2.3522365) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.464444) q[1];
sx q[1];
rz(-0.943999) q[1];
sx q[1];
rz(1.6897175) q[1];
rz(-3.0074189) q[3];
sx q[3];
rz(-0.65376702) q[3];
sx q[3];
rz(-1.8133481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0306586) q[2];
sx q[2];
rz(-1.3867644) q[2];
sx q[2];
rz(0.78732642) q[2];
rz(2.2287255) q[3];
sx q[3];
rz(-0.74936167) q[3];
sx q[3];
rz(-1.0095989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71516365) q[0];
sx q[0];
rz(-0.63868317) q[0];
sx q[0];
rz(-0.062967904) q[0];
rz(0.12403034) q[1];
sx q[1];
rz(-0.80563671) q[1];
sx q[1];
rz(-0.45809349) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8100909) q[0];
sx q[0];
rz(-1.1281663) q[0];
sx q[0];
rz(1.5682194) q[0];
x q[1];
rz(0.96180054) q[2];
sx q[2];
rz(-0.70448175) q[2];
sx q[2];
rz(1.5722164) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.304368) q[1];
sx q[1];
rz(-0.75290426) q[1];
sx q[1];
rz(1.7423082) q[1];
rz(-1.3985653) q[3];
sx q[3];
rz(-0.50695626) q[3];
sx q[3];
rz(0.34190049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1725585) q[2];
sx q[2];
rz(-2.2183552) q[2];
sx q[2];
rz(0.22949533) q[2];
rz(3.138792) q[3];
sx q[3];
rz(-2.2715748) q[3];
sx q[3];
rz(-1.3389448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5795508) q[0];
sx q[0];
rz(-2.8631449) q[0];
sx q[0];
rz(-0.91947412) q[0];
rz(-0.062285034) q[1];
sx q[1];
rz(-1.0039763) q[1];
sx q[1];
rz(1.2671635) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7125268) q[0];
sx q[0];
rz(-0.8493087) q[0];
sx q[0];
rz(2.6119786) q[0];
rz(0.15684776) q[2];
sx q[2];
rz(-0.75011293) q[2];
sx q[2];
rz(1.4175121) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8333203) q[1];
sx q[1];
rz(-0.38494021) q[1];
sx q[1];
rz(1.0233378) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1250161) q[3];
sx q[3];
rz(-2.6187754) q[3];
sx q[3];
rz(-1.0558053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1798114) q[2];
sx q[2];
rz(-0.87819019) q[2];
sx q[2];
rz(-0.58376694) q[2];
rz(-0.70872712) q[3];
sx q[3];
rz(-1.830359) q[3];
sx q[3];
rz(-3.1183929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0631183) q[0];
sx q[0];
rz(-0.45409504) q[0];
sx q[0];
rz(1.0725347) q[0];
rz(2.5947) q[1];
sx q[1];
rz(-1.8976338) q[1];
sx q[1];
rz(-2.0297208) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3277153) q[0];
sx q[0];
rz(-1.0881256) q[0];
sx q[0];
rz(3.0359603) q[0];
rz(0.56356168) q[2];
sx q[2];
rz(-2.5346018) q[2];
sx q[2];
rz(-0.8286455) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.48549451) q[1];
sx q[1];
rz(-1.494325) q[1];
sx q[1];
rz(2.8511726) q[1];
x q[2];
rz(2.7518919) q[3];
sx q[3];
rz(-1.1003564) q[3];
sx q[3];
rz(-2.5999992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.303858) q[2];
sx q[2];
rz(-1.4759109) q[2];
sx q[2];
rz(-1.1676577) q[2];
rz(1.5363103) q[3];
sx q[3];
rz(-1.6746018) q[3];
sx q[3];
rz(-1.7355828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4090356) q[0];
sx q[0];
rz(-1.5719825) q[0];
sx q[0];
rz(2.4107966) q[0];
rz(-0.90019512) q[1];
sx q[1];
rz(-2.3370445) q[1];
sx q[1];
rz(0.75497595) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2728111) q[0];
sx q[0];
rz(-2.4065354) q[0];
sx q[0];
rz(1.3520157) q[0];
rz(-2.412699) q[2];
sx q[2];
rz(-1.1155827) q[2];
sx q[2];
rz(-0.58065562) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4432115) q[1];
sx q[1];
rz(-1.6834007) q[1];
sx q[1];
rz(0.57971445) q[1];
rz(-pi) q[2];
rz(0.98324361) q[3];
sx q[3];
rz(-1.1508815) q[3];
sx q[3];
rz(-2.0411232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.76453152) q[2];
sx q[2];
rz(-1.7691282) q[2];
sx q[2];
rz(2.5047452) q[2];
rz(-2.8751255) q[3];
sx q[3];
rz(-1.0576495) q[3];
sx q[3];
rz(-1.5554957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72717845) q[0];
sx q[0];
rz(-1.1295015) q[0];
sx q[0];
rz(-2.0027347) q[0];
rz(0.75421929) q[1];
sx q[1];
rz(-0.33640877) q[1];
sx q[1];
rz(-3.1220904) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7050539) q[0];
sx q[0];
rz(-2.9199613) q[0];
sx q[0];
rz(-1.9428695) q[0];
rz(-pi) q[1];
rz(-0.15140622) q[2];
sx q[2];
rz(-1.7689686) q[2];
sx q[2];
rz(-2.1422447) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.49669493) q[1];
sx q[1];
rz(-0.95458191) q[1];
sx q[1];
rz(2.9195021) q[1];
rz(-pi) q[2];
rz(1.2583624) q[3];
sx q[3];
rz(-1.4133246) q[3];
sx q[3];
rz(-2.0563994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.0043682178) q[2];
sx q[2];
rz(-1.7251816) q[2];
sx q[2];
rz(-2.2926889) q[2];
rz(2.7539339) q[3];
sx q[3];
rz(-2.013423) q[3];
sx q[3];
rz(1.5415812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3432817) q[0];
sx q[0];
rz(-0.16769519) q[0];
sx q[0];
rz(2.6570901) q[0];
rz(-1.3867406) q[1];
sx q[1];
rz(-1.7157028) q[1];
sx q[1];
rz(-1.1482931) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.590811) q[0];
sx q[0];
rz(-1.5494487) q[0];
sx q[0];
rz(-0.14069964) q[0];
x q[1];
rz(-0.30015595) q[2];
sx q[2];
rz(-0.65320063) q[2];
sx q[2];
rz(2.4094827) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.82603589) q[1];
sx q[1];
rz(-2.450374) q[1];
sx q[1];
rz(-1.8408937) q[1];
rz(-pi) q[2];
rz(-0.34928068) q[3];
sx q[3];
rz(-1.619907) q[3];
sx q[3];
rz(1.1790566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.91223532) q[2];
sx q[2];
rz(-1.8476013) q[2];
sx q[2];
rz(-2.7764376) q[2];
rz(-0.12864104) q[3];
sx q[3];
rz(-1.235685) q[3];
sx q[3];
rz(-0.45583367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1098332) q[0];
sx q[0];
rz(-0.84072996) q[0];
sx q[0];
rz(1.6054556) q[0];
rz(0.96314349) q[1];
sx q[1];
rz(-1.8704725) q[1];
sx q[1];
rz(2.0830547) q[1];
rz(2.4776496) q[2];
sx q[2];
rz(-1.2524111) q[2];
sx q[2];
rz(-1.8128957) q[2];
rz(-0.036988463) q[3];
sx q[3];
rz(-2.130641) q[3];
sx q[3];
rz(-0.11161042) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
