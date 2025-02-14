OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2430159) q[0];
sx q[0];
rz(-1.1216811) q[0];
sx q[0];
rz(0.46410528) q[0];
rz(2.3857181) q[1];
sx q[1];
rz(2.2321489) q[1];
sx q[1];
rz(7.5997054) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1027604) q[0];
sx q[0];
rz(-1.505672) q[0];
sx q[0];
rz(-1.6054356) q[0];
rz(-1.6208956) q[2];
sx q[2];
rz(-1.9110137) q[2];
sx q[2];
rz(0.17776793) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.39151788) q[1];
sx q[1];
rz(-0.98528359) q[1];
sx q[1];
rz(1.2445009) q[1];
x q[2];
rz(2.2881704) q[3];
sx q[3];
rz(-2.7629921) q[3];
sx q[3];
rz(0.84747696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2509649) q[2];
sx q[2];
rz(-0.35718063) q[2];
sx q[2];
rz(2.8446021) q[2];
rz(2.4885528) q[3];
sx q[3];
rz(-1.5584471) q[3];
sx q[3];
rz(2.0853341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5121269) q[0];
sx q[0];
rz(-2.1194206) q[0];
sx q[0];
rz(-2.8566991) q[0];
rz(3.0902872) q[1];
sx q[1];
rz(-0.76301328) q[1];
sx q[1];
rz(-2.5047393) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013971) q[0];
sx q[0];
rz(-0.68249615) q[0];
sx q[0];
rz(2.3861814) q[0];
x q[1];
rz(3.1383855) q[2];
sx q[2];
rz(-0.35164552) q[2];
sx q[2];
rz(-0.52709792) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.53583685) q[1];
sx q[1];
rz(-2.9510289) q[1];
sx q[1];
rz(2.7298353) q[1];
rz(-pi) q[2];
rz(2.4356682) q[3];
sx q[3];
rz(-1.4691938) q[3];
sx q[3];
rz(1.3639579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.18179831) q[2];
sx q[2];
rz(-2.2836298) q[2];
sx q[2];
rz(0.96724969) q[2];
rz(-0.26425427) q[3];
sx q[3];
rz(-1.6382917) q[3];
sx q[3];
rz(1.2359469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0018175) q[0];
sx q[0];
rz(-1.5793261) q[0];
sx q[0];
rz(1.1199957) q[0];
rz(0.4758052) q[1];
sx q[1];
rz(-2.0640524) q[1];
sx q[1];
rz(2.8376104) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7426152) q[0];
sx q[0];
rz(-1.0053867) q[0];
sx q[0];
rz(3.1088367) q[0];
rz(-1.6725704) q[2];
sx q[2];
rz(-1.3341122) q[2];
sx q[2];
rz(2.9477811) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0416154) q[1];
sx q[1];
rz(-1.4681889) q[1];
sx q[1];
rz(2.2907738) q[1];
x q[2];
rz(-2.899053) q[3];
sx q[3];
rz(-2.4437196) q[3];
sx q[3];
rz(2.1062811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1659282) q[2];
sx q[2];
rz(-1.0378446) q[2];
sx q[2];
rz(0.69592875) q[2];
rz(-2.0078697) q[3];
sx q[3];
rz(-1.9418095) q[3];
sx q[3];
rz(-0.55293647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43941471) q[0];
sx q[0];
rz(-1.004383) q[0];
sx q[0];
rz(2.0942005) q[0];
rz(-1.9678496) q[1];
sx q[1];
rz(-0.67432299) q[1];
sx q[1];
rz(-1.5210927) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3701909) q[0];
sx q[0];
rz(-2.6242497) q[0];
sx q[0];
rz(2.8964554) q[0];
rz(-pi) q[1];
x q[1];
rz(0.32033605) q[2];
sx q[2];
rz(-0.89010677) q[2];
sx q[2];
rz(2.5706511) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2167705) q[1];
sx q[1];
rz(-0.7982648) q[1];
sx q[1];
rz(2.5096748) q[1];
rz(-pi) q[2];
rz(-1.8549881) q[3];
sx q[3];
rz(-1.8833369) q[3];
sx q[3];
rz(-0.72060637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0461222) q[2];
sx q[2];
rz(-1.0669402) q[2];
sx q[2];
rz(-2.6037237) q[2];
rz(-1.4872023) q[3];
sx q[3];
rz(-2.7345149) q[3];
sx q[3];
rz(0.68638221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9293514) q[0];
sx q[0];
rz(-0.23433267) q[0];
sx q[0];
rz(0.60337639) q[0];
rz(0.76250184) q[1];
sx q[1];
rz(-1.1926032) q[1];
sx q[1];
rz(0.71294436) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0423316) q[0];
sx q[0];
rz(-1.6368027) q[0];
sx q[0];
rz(-1.1724654) q[0];
rz(-1.2961757) q[2];
sx q[2];
rz(-2.6634187) q[2];
sx q[2];
rz(0.51152767) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0281369) q[1];
sx q[1];
rz(-1.077768) q[1];
sx q[1];
rz(2.9131469) q[1];
x q[2];
rz(-0.74728031) q[3];
sx q[3];
rz(-0.76939121) q[3];
sx q[3];
rz(-1.4667433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.9731628) q[2];
sx q[2];
rz(-0.91732401) q[2];
sx q[2];
rz(1.470835) q[2];
rz(2.6245608) q[3];
sx q[3];
rz(-1.0443338) q[3];
sx q[3];
rz(1.8487336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9789199) q[0];
sx q[0];
rz(-2.6277268) q[0];
sx q[0];
rz(2.5901219) q[0];
rz(-0.035471352) q[1];
sx q[1];
rz(-1.1330117) q[1];
sx q[1];
rz(1.5303401) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96824232) q[0];
sx q[0];
rz(-0.98778546) q[0];
sx q[0];
rz(0.75847404) q[0];
rz(-pi) q[1];
rz(1.6257203) q[2];
sx q[2];
rz(-0.74797219) q[2];
sx q[2];
rz(-2.2750281) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9236218) q[1];
sx q[1];
rz(-1.21232) q[1];
sx q[1];
rz(-2.3450407) q[1];
x q[2];
rz(2.2238268) q[3];
sx q[3];
rz(-1.6380596) q[3];
sx q[3];
rz(-1.3239087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0643206) q[2];
sx q[2];
rz(-0.52017009) q[2];
sx q[2];
rz(1.2426097) q[2];
rz(2.6284435) q[3];
sx q[3];
rz(-0.41392252) q[3];
sx q[3];
rz(-1.7665524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.649491) q[0];
sx q[0];
rz(-0.92925564) q[0];
sx q[0];
rz(-3.0249366) q[0];
rz(0.77313441) q[1];
sx q[1];
rz(-1.1583068) q[1];
sx q[1];
rz(-1.6366417) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9336108) q[0];
sx q[0];
rz(-0.79219063) q[0];
sx q[0];
rz(-1.1536145) q[0];
rz(-pi) q[1];
rz(3.0819986) q[2];
sx q[2];
rz(-1.182654) q[2];
sx q[2];
rz(1.3898894) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.22681397) q[1];
sx q[1];
rz(-0.89229167) q[1];
sx q[1];
rz(-2.6752276) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2833474) q[3];
sx q[3];
rz(-0.93837591) q[3];
sx q[3];
rz(2.8638864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.92631212) q[2];
sx q[2];
rz(-2.8947688) q[2];
sx q[2];
rz(-0.8872633) q[2];
rz(0.67982802) q[3];
sx q[3];
rz(-0.66893783) q[3];
sx q[3];
rz(0.63290709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.585084) q[0];
sx q[0];
rz(-1.8483138) q[0];
sx q[0];
rz(-2.4853117) q[0];
rz(2.9346924) q[1];
sx q[1];
rz(-1.2196536) q[1];
sx q[1];
rz(1.9237178) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6164056) q[0];
sx q[0];
rz(-1.5357247) q[0];
sx q[0];
rz(-0.082790815) q[0];
rz(-pi) q[1];
rz(-0.98274173) q[2];
sx q[2];
rz(-1.2763192) q[2];
sx q[2];
rz(-0.5628995) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.52384842) q[1];
sx q[1];
rz(-2.9459125) q[1];
sx q[1];
rz(-2.5168672) q[1];
x q[2];
rz(0.37533203) q[3];
sx q[3];
rz(-2.4960244) q[3];
sx q[3];
rz(0.69824346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3394341) q[2];
sx q[2];
rz(-1.3863486) q[2];
sx q[2];
rz(0.68186861) q[2];
rz(1.1784461) q[3];
sx q[3];
rz(-0.75740564) q[3];
sx q[3];
rz(1.9044378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39663974) q[0];
sx q[0];
rz(-1.5322026) q[0];
sx q[0];
rz(-2.9314281) q[0];
rz(0.22178966) q[1];
sx q[1];
rz(-2.2237325) q[1];
sx q[1];
rz(-0.04235696) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1975178) q[0];
sx q[0];
rz(-2.952708) q[0];
sx q[0];
rz(-3.023639) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3826959) q[2];
sx q[2];
rz(-0.79115552) q[2];
sx q[2];
rz(-3.0393785) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.48407979) q[1];
sx q[1];
rz(-1.8575605) q[1];
sx q[1];
rz(-0.57876719) q[1];
rz(0.5333535) q[3];
sx q[3];
rz(-1.9755529) q[3];
sx q[3];
rz(-3.0572756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.28045851) q[2];
sx q[2];
rz(-2.0529842) q[2];
sx q[2];
rz(-2.4750278) q[2];
rz(2.376453) q[3];
sx q[3];
rz(-2.9355526) q[3];
sx q[3];
rz(1.4008745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6078981) q[0];
sx q[0];
rz(-0.38132897) q[0];
sx q[0];
rz(-2.4545942) q[0];
rz(-1.2835361) q[1];
sx q[1];
rz(-1.9891519) q[1];
sx q[1];
rz(-2.5680465) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4939369) q[0];
sx q[0];
rz(-2.082061) q[0];
sx q[0];
rz(2.4535116) q[0];
rz(2.9098689) q[2];
sx q[2];
rz(-2.0123693) q[2];
sx q[2];
rz(-1.9927466) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6489689) q[1];
sx q[1];
rz(-2.0855806) q[1];
sx q[1];
rz(2.749498) q[1];
x q[2];
rz(-2.2881094) q[3];
sx q[3];
rz(-2.1721345) q[3];
sx q[3];
rz(-0.36313148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0428697) q[2];
sx q[2];
rz(-1.9645773) q[2];
sx q[2];
rz(-3.0549808) q[2];
rz(-3.0232271) q[3];
sx q[3];
rz(-1.5990853) q[3];
sx q[3];
rz(0.043244403) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31556986) q[0];
sx q[0];
rz(-2.12349) q[0];
sx q[0];
rz(0.77793599) q[0];
rz(2.8553873) q[1];
sx q[1];
rz(-2.310391) q[1];
sx q[1];
rz(-1.8117767) q[1];
rz(2.8092842) q[2];
sx q[2];
rz(-1.7996017) q[2];
sx q[2];
rz(-1.9881291) q[2];
rz(0.44888857) q[3];
sx q[3];
rz(-1.8406624) q[3];
sx q[3];
rz(0.016907666) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
