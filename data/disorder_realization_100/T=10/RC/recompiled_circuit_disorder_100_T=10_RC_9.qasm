OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.93958062) q[0];
sx q[0];
rz(2.7913845) q[0];
sx q[0];
rz(12.933001) q[0];
rz(0.86759138) q[1];
sx q[1];
rz(-2.4974513) q[1];
sx q[1];
rz(-1.6860513) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7830084) q[0];
sx q[0];
rz(-1.2281679) q[0];
sx q[0];
rz(1.1230099) q[0];
rz(-pi) q[1];
rz(0.86906616) q[2];
sx q[2];
rz(-1.7742187) q[2];
sx q[2];
rz(-2.2609401) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3914324) q[1];
sx q[1];
rz(-0.62392985) q[1];
sx q[1];
rz(-2.0425914) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0936712) q[3];
sx q[3];
rz(-1.5600292) q[3];
sx q[3];
rz(2.3445232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.036844) q[2];
sx q[2];
rz(-1.8061183) q[2];
sx q[2];
rz(2.74995) q[2];
rz(3.0018905) q[3];
sx q[3];
rz(-0.67290664) q[3];
sx q[3];
rz(3.1203549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4166819) q[0];
sx q[0];
rz(-1.7371438) q[0];
sx q[0];
rz(1.8649944) q[0];
rz(2.2712767) q[1];
sx q[1];
rz(-1.566193) q[1];
sx q[1];
rz(-1.93719) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8613974) q[0];
sx q[0];
rz(-1.9635927) q[0];
sx q[0];
rz(2.0307721) q[0];
rz(-pi) q[1];
rz(-2.9559829) q[2];
sx q[2];
rz(-2.2679272) q[2];
sx q[2];
rz(-2.6999161) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7960297) q[1];
sx q[1];
rz(-1.2064762) q[1];
sx q[1];
rz(-1.6771392) q[1];
x q[2];
rz(0.28537206) q[3];
sx q[3];
rz(-2.223569) q[3];
sx q[3];
rz(-0.6245581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.63699547) q[2];
sx q[2];
rz(-2.6837139) q[2];
sx q[2];
rz(1.9223928) q[2];
rz(0.35456625) q[3];
sx q[3];
rz(-0.48186007) q[3];
sx q[3];
rz(-1.5725117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59042674) q[0];
sx q[0];
rz(-1.6141491) q[0];
sx q[0];
rz(-2.5235126) q[0];
rz(-0.016050054) q[1];
sx q[1];
rz(-2.2647104) q[1];
sx q[1];
rz(-1.9504257) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58142692) q[0];
sx q[0];
rz(-1.418857) q[0];
sx q[0];
rz(-1.2537969) q[0];
rz(0.70706681) q[2];
sx q[2];
rz(-0.36598772) q[2];
sx q[2];
rz(-2.6235839) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1515854) q[1];
sx q[1];
rz(-1.7813492) q[1];
sx q[1];
rz(-1.7819808) q[1];
x q[2];
rz(-2.5268764) q[3];
sx q[3];
rz(-2.5066272) q[3];
sx q[3];
rz(0.8599417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.96737343) q[2];
sx q[2];
rz(-2.1790395) q[2];
sx q[2];
rz(-2.1276786) q[2];
rz(-1.3570471) q[3];
sx q[3];
rz(-2.3116528) q[3];
sx q[3];
rz(-1.0323662) q[3];
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
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8106666) q[0];
sx q[0];
rz(-1.4241011) q[0];
sx q[0];
rz(-0.20430918) q[0];
rz(-1.7640242) q[1];
sx q[1];
rz(-1.7267449) q[1];
sx q[1];
rz(-0.92299443) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8078634) q[0];
sx q[0];
rz(-2.6834052) q[0];
sx q[0];
rz(0.74497594) q[0];
x q[1];
rz(1.5590018) q[2];
sx q[2];
rz(-1.5268541) q[2];
sx q[2];
rz(-2.5269788) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.97913617) q[1];
sx q[1];
rz(-0.41185954) q[1];
sx q[1];
rz(-1.3084175) q[1];
x q[2];
rz(3.119167) q[3];
sx q[3];
rz(-2.2750345) q[3];
sx q[3];
rz(-1.1426329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6874281) q[2];
sx q[2];
rz(-1.5565846) q[2];
sx q[2];
rz(-0.50951177) q[2];
rz(2.4984958) q[3];
sx q[3];
rz(-2.0431079) q[3];
sx q[3];
rz(2.174214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61280695) q[0];
sx q[0];
rz(-1.9354154) q[0];
sx q[0];
rz(0.91941419) q[0];
rz(0.98584229) q[1];
sx q[1];
rz(-1.3580094) q[1];
sx q[1];
rz(1.3607508) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1795792) q[0];
sx q[0];
rz(-2.8498581) q[0];
sx q[0];
rz(-2.4289262) q[0];
rz(-pi) q[1];
rz(-2.1941357) q[2];
sx q[2];
rz(-0.91295487) q[2];
sx q[2];
rz(2.5156977) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.67140019) q[1];
sx q[1];
rz(-0.50506401) q[1];
sx q[1];
rz(-2.7027674) q[1];
rz(-pi) q[2];
rz(-2.1702483) q[3];
sx q[3];
rz(-1.6347307) q[3];
sx q[3];
rz(-1.5980699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3551066) q[2];
sx q[2];
rz(-1.7559218) q[2];
sx q[2];
rz(-0.11486593) q[2];
rz(-0.45587513) q[3];
sx q[3];
rz(-1.3331648) q[3];
sx q[3];
rz(1.8615287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24546656) q[0];
sx q[0];
rz(-1.2446612) q[0];
sx q[0];
rz(-1.2448357) q[0];
rz(2.4339829) q[1];
sx q[1];
rz(-1.6376303) q[1];
sx q[1];
rz(0.27522603) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48588612) q[0];
sx q[0];
rz(-2.5285892) q[0];
sx q[0];
rz(-0.97047752) q[0];
rz(-0.94324865) q[2];
sx q[2];
rz(-2.3970282) q[2];
sx q[2];
rz(-2.9668691) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3616398) q[1];
sx q[1];
rz(-2.0556903) q[1];
sx q[1];
rz(1.9655025) q[1];
rz(-pi) q[2];
rz(-0.082645881) q[3];
sx q[3];
rz(-2.189889) q[3];
sx q[3];
rz(-1.1766528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.747921) q[2];
sx q[2];
rz(-1.5449521) q[2];
sx q[2];
rz(0.45864027) q[2];
rz(-0.7115055) q[3];
sx q[3];
rz(-0.68370521) q[3];
sx q[3];
rz(2.5512364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28073072) q[0];
sx q[0];
rz(-1.9545398) q[0];
sx q[0];
rz(1.77805) q[0];
rz(1.7215615) q[1];
sx q[1];
rz(-2.6224711) q[1];
sx q[1];
rz(-2.4218959) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62455432) q[0];
sx q[0];
rz(-1.2861684) q[0];
sx q[0];
rz(0.65529234) q[0];
x q[1];
rz(2.2816706) q[2];
sx q[2];
rz(-0.93647525) q[2];
sx q[2];
rz(0.89564043) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3230348) q[1];
sx q[1];
rz(-2.5948988) q[1];
sx q[1];
rz(-0.085303765) q[1];
x q[2];
rz(-0.32960524) q[3];
sx q[3];
rz(-0.70809396) q[3];
sx q[3];
rz(3.0081188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.482243) q[2];
sx q[2];
rz(-1.6094001) q[2];
sx q[2];
rz(-2.4534295) q[2];
rz(0.041953772) q[3];
sx q[3];
rz(-2.0418906) q[3];
sx q[3];
rz(0.28013128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91350895) q[0];
sx q[0];
rz(-1.9788454) q[0];
sx q[0];
rz(2.9550609) q[0];
rz(-0.60449156) q[1];
sx q[1];
rz(-1.0124606) q[1];
sx q[1];
rz(1.4950745) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.625531) q[0];
sx q[0];
rz(-0.4058668) q[0];
sx q[0];
rz(-2.5748475) q[0];
rz(-pi) q[1];
x q[1];
rz(1.375884) q[2];
sx q[2];
rz(-0.88390985) q[2];
sx q[2];
rz(1.9538823) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9547792) q[1];
sx q[1];
rz(-0.28007945) q[1];
sx q[1];
rz(2.7571452) q[1];
x q[2];
rz(2.5008194) q[3];
sx q[3];
rz(-2.2230004) q[3];
sx q[3];
rz(0.81354248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.24215332) q[2];
sx q[2];
rz(-2.184325) q[2];
sx q[2];
rz(-3.126826) q[2];
rz(1.2773369) q[3];
sx q[3];
rz(-1.3093964) q[3];
sx q[3];
rz(1.7285255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16167851) q[0];
sx q[0];
rz(-2.4665687) q[0];
sx q[0];
rz(2.263608) q[0];
rz(2.9453078) q[1];
sx q[1];
rz(-1.1801964) q[1];
sx q[1];
rz(1.3605798) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6920647) q[0];
sx q[0];
rz(-0.99943107) q[0];
sx q[0];
rz(-3.0968303) q[0];
rz(0.94366818) q[2];
sx q[2];
rz(-1.5329648) q[2];
sx q[2];
rz(-0.19247069) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9578743) q[1];
sx q[1];
rz(-2.6596333) q[1];
sx q[1];
rz(-2.3493489) q[1];
x q[2];
rz(-1.1514879) q[3];
sx q[3];
rz(-1.7125704) q[3];
sx q[3];
rz(-1.7820953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.517841) q[2];
sx q[2];
rz(-1.6343642) q[2];
sx q[2];
rz(-2.2488135) q[2];
rz(-3.1023074) q[3];
sx q[3];
rz(-1.4792484) q[3];
sx q[3];
rz(-2.3915496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84353012) q[0];
sx q[0];
rz(-0.003304464) q[0];
sx q[0];
rz(3.0950586) q[0];
rz(2.2325366) q[1];
sx q[1];
rz(-0.87288705) q[1];
sx q[1];
rz(0.7199026) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7881308) q[0];
sx q[0];
rz(-0.25665584) q[0];
sx q[0];
rz(-1.8887595) q[0];
rz(1.5653531) q[2];
sx q[2];
rz(-1.2619702) q[2];
sx q[2];
rz(0.18143166) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9626999) q[1];
sx q[1];
rz(-1.8176515) q[1];
sx q[1];
rz(3.1256413) q[1];
x q[2];
rz(2.3941444) q[3];
sx q[3];
rz(-0.8851074) q[3];
sx q[3];
rz(1.7789343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4204734) q[2];
sx q[2];
rz(-1.3995918) q[2];
sx q[2];
rz(0.026467888) q[2];
rz(1.8376393) q[3];
sx q[3];
rz(-0.81726685) q[3];
sx q[3];
rz(1.2560237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7883041) q[0];
sx q[0];
rz(-1.7699387) q[0];
sx q[0];
rz(1.8040245) q[0];
rz(-2.9251255) q[1];
sx q[1];
rz(-1.4420061) q[1];
sx q[1];
rz(-1.9056086) q[1];
rz(-2.4344865) q[2];
sx q[2];
rz(-1.3277935) q[2];
sx q[2];
rz(-2.1869616) q[2];
rz(1.5218232) q[3];
sx q[3];
rz(-2.5544142) q[3];
sx q[3];
rz(2.1094473) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];