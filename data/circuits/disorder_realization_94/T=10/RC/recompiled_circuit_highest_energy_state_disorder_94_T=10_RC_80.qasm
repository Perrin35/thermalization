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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038832207) q[0];
sx q[0];
rz(-1.6359207) q[0];
sx q[0];
rz(-1.536157) q[0];
rz(-pi) q[1];
rz(2.8009802) q[2];
sx q[2];
rz(-1.5235708) q[2];
sx q[2];
rz(1.7318341) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1471828) q[1];
sx q[1];
rz(-1.8412245) q[1];
sx q[1];
rz(2.5308591) q[1];
rz(-pi) q[2];
rz(-0.85342225) q[3];
sx q[3];
rz(-0.3786006) q[3];
sx q[3];
rz(-0.84747696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2509649) q[2];
sx q[2];
rz(-0.35718063) q[2];
sx q[2];
rz(2.8446021) q[2];
rz(-2.4885528) q[3];
sx q[3];
rz(-1.5831455) q[3];
sx q[3];
rz(-1.0562586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5121269) q[0];
sx q[0];
rz(-2.1194206) q[0];
sx q[0];
rz(2.8566991) q[0];
rz(-0.051305436) q[1];
sx q[1];
rz(-2.3785794) q[1];
sx q[1];
rz(2.5047393) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2460265) q[0];
sx q[0];
rz(-2.0478529) q[0];
sx q[0];
rz(-2.0791847) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0032071438) q[2];
sx q[2];
rz(-0.35164552) q[2];
sx q[2];
rz(2.6144947) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7015345) q[1];
sx q[1];
rz(-1.6466759) q[1];
sx q[1];
rz(-2.9666191) q[1];
rz(0.15588197) q[3];
sx q[3];
rz(-2.4296399) q[3];
sx q[3];
rz(-0.32526325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9597943) q[2];
sx q[2];
rz(-0.85796285) q[2];
sx q[2];
rz(-0.96724969) q[2];
rz(-0.26425427) q[3];
sx q[3];
rz(-1.5033009) q[3];
sx q[3];
rz(1.9056457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13977519) q[0];
sx q[0];
rz(-1.5793261) q[0];
sx q[0];
rz(-1.1199957) q[0];
rz(-0.4758052) q[1];
sx q[1];
rz(-2.0640524) q[1];
sx q[1];
rz(0.30398223) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39897746) q[0];
sx q[0];
rz(-1.0053867) q[0];
sx q[0];
rz(3.1088367) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4690222) q[2];
sx q[2];
rz(-1.3341122) q[2];
sx q[2];
rz(2.9477811) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5811823) q[1];
sx q[1];
rz(-0.85542233) q[1];
sx q[1];
rz(-0.13611273) q[1];
x q[2];
rz(-2.4583057) q[3];
sx q[3];
rz(-1.7257455) q[3];
sx q[3];
rz(-2.4187627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9756644) q[2];
sx q[2];
rz(-2.1037481) q[2];
sx q[2];
rz(2.4456639) q[2];
rz(1.1337229) q[3];
sx q[3];
rz(-1.1997831) q[3];
sx q[3];
rz(0.55293647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7021779) q[0];
sx q[0];
rz(-2.1372097) q[0];
sx q[0];
rz(-2.0942005) q[0];
rz(1.1737431) q[1];
sx q[1];
rz(-0.67432299) q[1];
sx q[1];
rz(-1.5210927) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58529638) q[0];
sx q[0];
rz(-1.4504787) q[0];
sx q[0];
rz(2.6371945) q[0];
x q[1];
rz(0.32033605) q[2];
sx q[2];
rz(-2.2514859) q[2];
sx q[2];
rz(-2.5706511) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7340659) q[1];
sx q[1];
rz(-0.95469771) q[1];
sx q[1];
rz(1.0259088) q[1];
rz(-pi) q[2];
x q[2];
rz(0.71468227) q[3];
sx q[3];
rz(-2.7223248) q[3];
sx q[3];
rz(-1.4803606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.095470458) q[2];
sx q[2];
rz(-1.0669402) q[2];
sx q[2];
rz(0.53786892) q[2];
rz(1.6543903) q[3];
sx q[3];
rz(-0.40707773) q[3];
sx q[3];
rz(-0.68638221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9293514) q[0];
sx q[0];
rz(-2.90726) q[0];
sx q[0];
rz(2.5382163) q[0];
rz(2.3790908) q[1];
sx q[1];
rz(-1.9489894) q[1];
sx q[1];
rz(-2.4286483) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099261053) q[0];
sx q[0];
rz(-1.5047899) q[0];
sx q[0];
rz(1.9691273) q[0];
rz(-pi) q[1];
rz(1.845417) q[2];
sx q[2];
rz(-0.47817395) q[2];
sx q[2];
rz(2.630065) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5746497) q[1];
sx q[1];
rz(-1.7716367) q[1];
sx q[1];
rz(-1.0667136) q[1];
rz(-0.98865786) q[3];
sx q[3];
rz(-2.1063559) q[3];
sx q[3];
rz(0.55544191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.9731628) q[2];
sx q[2];
rz(-0.91732401) q[2];
sx q[2];
rz(1.6707576) q[2];
rz(-2.6245608) q[3];
sx q[3];
rz(-2.0972589) q[3];
sx q[3];
rz(-1.2928591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1626728) q[0];
sx q[0];
rz(-2.6277268) q[0];
sx q[0];
rz(2.5901219) q[0];
rz(0.035471352) q[1];
sx q[1];
rz(-1.1330117) q[1];
sx q[1];
rz(1.6112526) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0198675) q[0];
sx q[0];
rz(-2.1824153) q[0];
sx q[0];
rz(-2.3082971) q[0];
rz(-pi) q[1];
rz(-3.0907029) q[2];
sx q[2];
rz(-2.3173703) q[2];
sx q[2];
rz(0.79170601) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9236218) q[1];
sx q[1];
rz(-1.9292726) q[1];
sx q[1];
rz(0.79655194) q[1];
rz(-pi) q[2];
rz(0.91776589) q[3];
sx q[3];
rz(-1.5035331) q[3];
sx q[3];
rz(1.8176839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0772721) q[2];
sx q[2];
rz(-0.52017009) q[2];
sx q[2];
rz(1.8989829) q[2];
rz(2.6284435) q[3];
sx q[3];
rz(-0.41392252) q[3];
sx q[3];
rz(1.3750403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.649491) q[0];
sx q[0];
rz(-0.92925564) q[0];
sx q[0];
rz(0.11665601) q[0];
rz(-0.77313441) q[1];
sx q[1];
rz(-1.9832858) q[1];
sx q[1];
rz(1.504951) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0805252) q[0];
sx q[0];
rz(-1.2781898) q[0];
sx q[0];
rz(-0.82340296) q[0];
rz(-1.9595615) q[2];
sx q[2];
rz(-1.6259527) q[2];
sx q[2];
rz(2.9832632) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.44918647) q[1];
sx q[1];
rz(-2.3397361) q[1];
sx q[1];
rz(-1.0620326) q[1];
rz(-pi) q[2];
rz(-1.2833474) q[3];
sx q[3];
rz(-2.2032167) q[3];
sx q[3];
rz(-2.8638864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2152805) q[2];
sx q[2];
rz(-2.8947688) q[2];
sx q[2];
rz(2.2543294) q[2];
rz(2.4617646) q[3];
sx q[3];
rz(-2.4726548) q[3];
sx q[3];
rz(-2.5086856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.585084) q[0];
sx q[0];
rz(-1.2932788) q[0];
sx q[0];
rz(2.4853117) q[0];
rz(-2.9346924) q[1];
sx q[1];
rz(-1.2196536) q[1];
sx q[1];
rz(-1.9237178) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.098893) q[0];
sx q[0];
rz(-1.6535361) q[0];
sx q[0];
rz(-1.6059884) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1588509) q[2];
sx q[2];
rz(-1.2763192) q[2];
sx q[2];
rz(2.5786932) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9838502) q[1];
sx q[1];
rz(-1.4124253) q[1];
sx q[1];
rz(-1.4553797) q[1];
rz(-pi) q[2];
rz(-1.3013873) q[3];
sx q[3];
rz(-2.1649033) q[3];
sx q[3];
rz(1.9850933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8021585) q[2];
sx q[2];
rz(-1.3863486) q[2];
sx q[2];
rz(-2.459724) q[2];
rz(1.9631466) q[3];
sx q[3];
rz(-0.75740564) q[3];
sx q[3];
rz(-1.9044378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7449529) q[0];
sx q[0];
rz(-1.5322026) q[0];
sx q[0];
rz(-0.21016453) q[0];
rz(-0.22178966) q[1];
sx q[1];
rz(-2.2237325) q[1];
sx q[1];
rz(-3.0992357) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0641441) q[0];
sx q[0];
rz(-1.7583529) q[0];
sx q[0];
rz(-1.5932887) q[0];
x q[1];
rz(2.508411) q[2];
sx q[2];
rz(-2.0821619) q[2];
sx q[2];
rz(0.88054576) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4952462) q[1];
sx q[1];
rz(-0.63856417) q[1];
sx q[1];
rz(-2.6471443) q[1];
x q[2];
rz(-0.5333535) q[3];
sx q[3];
rz(-1.9755529) q[3];
sx q[3];
rz(3.0572756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.28045851) q[2];
sx q[2];
rz(-2.0529842) q[2];
sx q[2];
rz(-0.66656485) q[2];
rz(0.76513964) q[3];
sx q[3];
rz(-2.9355526) q[3];
sx q[3];
rz(-1.4008745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6078981) q[0];
sx q[0];
rz(-0.38132897) q[0];
sx q[0];
rz(0.68699849) q[0];
rz(1.2835361) q[1];
sx q[1];
rz(-1.9891519) q[1];
sx q[1];
rz(-0.57354617) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6015778) q[0];
sx q[0];
rz(-2.3099843) q[0];
sx q[0];
rz(0.7235781) q[0];
rz(-pi) q[1];
rz(1.1185455) q[2];
sx q[2];
rz(-0.49511038) q[2];
sx q[2];
rz(0.64436382) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0189218) q[1];
sx q[1];
rz(-1.9098567) q[1];
sx q[1];
rz(1.0215205) q[1];
rz(0.76400842) q[3];
sx q[3];
rz(-0.90029085) q[3];
sx q[3];
rz(-1.7830199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0428697) q[2];
sx q[2];
rz(-1.1770153) q[2];
sx q[2];
rz(-0.086611835) q[2];
rz(3.0232271) q[3];
sx q[3];
rz(-1.5990853) q[3];
sx q[3];
rz(3.0983483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31556986) q[0];
sx q[0];
rz(-1.0181027) q[0];
sx q[0];
rz(-2.3636567) q[0];
rz(-2.8553873) q[1];
sx q[1];
rz(-0.83120167) q[1];
sx q[1];
rz(1.3298159) q[1];
rz(1.329245) q[2];
sx q[2];
rz(-1.2474682) q[2];
sx q[2];
rz(-0.4954485) q[2];
rz(2.574118) q[3];
sx q[3];
rz(-0.51898659) q[3];
sx q[3];
rz(1.082194) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
