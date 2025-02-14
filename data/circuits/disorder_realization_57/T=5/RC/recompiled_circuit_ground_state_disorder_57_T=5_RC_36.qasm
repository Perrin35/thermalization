OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.39855555) q[0];
sx q[0];
rz(-0.86657137) q[0];
sx q[0];
rz(0.33696365) q[0];
rz(0.26951867) q[1];
sx q[1];
rz(4.0842847) q[1];
sx q[1];
rz(10.684291) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.082734) q[0];
sx q[0];
rz(-1.0089416) q[0];
sx q[0];
rz(0.64942067) q[0];
rz(1.7808229) q[2];
sx q[2];
rz(-2.2968596) q[2];
sx q[2];
rz(-3.0371916) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9971157) q[1];
sx q[1];
rz(-1.0273178) q[1];
sx q[1];
rz(-1.687084) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1416177) q[3];
sx q[3];
rz(-0.68203841) q[3];
sx q[3];
rz(1.2004971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.17244615) q[2];
sx q[2];
rz(-1.5643876) q[2];
sx q[2];
rz(-1.2338314) q[2];
rz(1.5305758) q[3];
sx q[3];
rz(-1.4065892) q[3];
sx q[3];
rz(2.5764901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3259786) q[0];
sx q[0];
rz(-2.1765206) q[0];
sx q[0];
rz(-2.8976029) q[0];
rz(1.1583534) q[1];
sx q[1];
rz(-1.0141677) q[1];
sx q[1];
rz(-0.44833952) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8566003) q[0];
sx q[0];
rz(-2.0584848) q[0];
sx q[0];
rz(1.3707121) q[0];
rz(-pi) q[1];
rz(2.3116287) q[2];
sx q[2];
rz(-1.1051851) q[2];
sx q[2];
rz(0.59620406) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.068591) q[1];
sx q[1];
rz(-1.0319971) q[1];
sx q[1];
rz(-0.52937845) q[1];
x q[2];
rz(-1.0226986) q[3];
sx q[3];
rz(-1.1493936) q[3];
sx q[3];
rz(2.2157026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1925194) q[2];
sx q[2];
rz(-0.039704617) q[2];
sx q[2];
rz(-0.81083361) q[2];
rz(0.0090946322) q[3];
sx q[3];
rz(-1.6154715) q[3];
sx q[3];
rz(-2.8240589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4119165) q[0];
sx q[0];
rz(-1.7971973) q[0];
sx q[0];
rz(2.1938531) q[0];
rz(-2.2432227) q[1];
sx q[1];
rz(-2.4285474) q[1];
sx q[1];
rz(2.9352442) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2437993) q[0];
sx q[0];
rz(-0.81598213) q[0];
sx q[0];
rz(-3.131098) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3487665) q[2];
sx q[2];
rz(-1.2812876) q[2];
sx q[2];
rz(-0.0057366554) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.84103407) q[1];
sx q[1];
rz(-1.7766469) q[1];
sx q[1];
rz(1.1080075) q[1];
x q[2];
rz(-0.32636362) q[3];
sx q[3];
rz(-1.6152528) q[3];
sx q[3];
rz(1.8281368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.72948939) q[2];
sx q[2];
rz(-0.99200839) q[2];
sx q[2];
rz(-0.98423973) q[2];
rz(2.6442243) q[3];
sx q[3];
rz(-1.8822742) q[3];
sx q[3];
rz(-2.3965912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7797778) q[0];
sx q[0];
rz(-0.091876939) q[0];
sx q[0];
rz(0.20633695) q[0];
rz(-1.4641209) q[1];
sx q[1];
rz(-0.76281491) q[1];
sx q[1];
rz(0.62613553) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7595644) q[0];
sx q[0];
rz(-0.81597486) q[0];
sx q[0];
rz(-1.0797281) q[0];
rz(0.10278662) q[2];
sx q[2];
rz(-1.572552) q[2];
sx q[2];
rz(2.2660915) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4724222) q[1];
sx q[1];
rz(-1.0423933) q[1];
sx q[1];
rz(1.0086568) q[1];
x q[2];
rz(0.098931626) q[3];
sx q[3];
rz(-1.9557992) q[3];
sx q[3];
rz(1.0855293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.76554406) q[2];
sx q[2];
rz(-0.77771336) q[2];
sx q[2];
rz(0.98975873) q[2];
rz(1.0342213) q[3];
sx q[3];
rz(-1.1690305) q[3];
sx q[3];
rz(-0.52556747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.702221) q[0];
sx q[0];
rz(-1.419743) q[0];
sx q[0];
rz(1.9787582) q[0];
rz(2.974466) q[1];
sx q[1];
rz(-1.6163328) q[1];
sx q[1];
rz(-2.4662245) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4136001) q[0];
sx q[0];
rz(-2.4420847) q[0];
sx q[0];
rz(-2.7986155) q[0];
rz(-pi) q[1];
rz(2.295601) q[2];
sx q[2];
rz(-0.48946871) q[2];
sx q[2];
rz(-2.6014525) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.45686238) q[1];
sx q[1];
rz(-0.56512512) q[1];
sx q[1];
rz(-0.88688382) q[1];
x q[2];
rz(0.97978239) q[3];
sx q[3];
rz(-1.3434935) q[3];
sx q[3];
rz(1.0349555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2165788) q[2];
sx q[2];
rz(-1.4224956) q[2];
sx q[2];
rz(-0.50755802) q[2];
rz(-3.1223068) q[3];
sx q[3];
rz(-2.3465893) q[3];
sx q[3];
rz(-1.9222586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5386706) q[0];
sx q[0];
rz(-0.62748533) q[0];
sx q[0];
rz(-0.70166171) q[0];
rz(2.498846) q[1];
sx q[1];
rz(-1.1593436) q[1];
sx q[1];
rz(1.3210375) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3468332) q[0];
sx q[0];
rz(-2.1982895) q[0];
sx q[0];
rz(-0.015608257) q[0];
x q[1];
rz(0.3568383) q[2];
sx q[2];
rz(-0.90274397) q[2];
sx q[2];
rz(2.2989863) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5112662) q[1];
sx q[1];
rz(-2.0980586) q[1];
sx q[1];
rz(0.43095891) q[1];
x q[2];
rz(-1.2679638) q[3];
sx q[3];
rz(-2.3583671) q[3];
sx q[3];
rz(2.2719605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.295149) q[2];
sx q[2];
rz(-0.48131338) q[2];
sx q[2];
rz(0.63977891) q[2];
rz(2.4844737) q[3];
sx q[3];
rz(-1.492447) q[3];
sx q[3];
rz(0.28321987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.039728634) q[0];
sx q[0];
rz(-1.855408) q[0];
sx q[0];
rz(2.3611327) q[0];
rz(-1.2377493) q[1];
sx q[1];
rz(-1.2982118) q[1];
sx q[1];
rz(2.9775528) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.964351) q[0];
sx q[0];
rz(-0.57946396) q[0];
sx q[0];
rz(2.01418) q[0];
rz(-pi) q[1];
rz(-0.82752821) q[2];
sx q[2];
rz(-2.568576) q[2];
sx q[2];
rz(2.1967251) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.32719993) q[1];
sx q[1];
rz(-0.76451028) q[1];
sx q[1];
rz(-0.79245497) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6884671) q[3];
sx q[3];
rz(-2.2966566) q[3];
sx q[3];
rz(0.22900861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.26874545) q[2];
sx q[2];
rz(-1.9445323) q[2];
sx q[2];
rz(-0.8806814) q[2];
rz(1.4745447) q[3];
sx q[3];
rz(-1.0677974) q[3];
sx q[3];
rz(-1.7269945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47377652) q[0];
sx q[0];
rz(-0.18874636) q[0];
sx q[0];
rz(-2.012398) q[0];
rz(-2.1881564) q[1];
sx q[1];
rz(-0.9402746) q[1];
sx q[1];
rz(2.5530691) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9669707) q[0];
sx q[0];
rz(-1.0292639) q[0];
sx q[0];
rz(-0.35090272) q[0];
rz(-1.2382201) q[2];
sx q[2];
rz(-1.3540478) q[2];
sx q[2];
rz(0.34453604) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5626763) q[1];
sx q[1];
rz(-2.2136218) q[1];
sx q[1];
rz(-1.972727) q[1];
rz(-pi) q[2];
rz(-3.088589) q[3];
sx q[3];
rz(-1.7946968) q[3];
sx q[3];
rz(2.8877024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9863161) q[2];
sx q[2];
rz(-1.0078112) q[2];
sx q[2];
rz(-2.0493832) q[2];
rz(-2.3327667) q[3];
sx q[3];
rz(-1.0816962) q[3];
sx q[3];
rz(1.697418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5876193) q[0];
sx q[0];
rz(-2.0646136) q[0];
sx q[0];
rz(0.67163604) q[0];
rz(-2.6486168) q[1];
sx q[1];
rz(-0.88075811) q[1];
sx q[1];
rz(-1.902045) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.400378) q[0];
sx q[0];
rz(-1.5733871) q[0];
sx q[0];
rz(1.8080312) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9443594) q[2];
sx q[2];
rz(-1.5805564) q[2];
sx q[2];
rz(-0.41744864) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0555435) q[1];
sx q[1];
rz(-0.99850875) q[1];
sx q[1];
rz(1.0466762) q[1];
rz(1.4080234) q[3];
sx q[3];
rz(-2.3331169) q[3];
sx q[3];
rz(-2.2089434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0121743) q[2];
sx q[2];
rz(-1.6678026) q[2];
sx q[2];
rz(2.5950281) q[2];
rz(-0.92285815) q[3];
sx q[3];
rz(-0.62896362) q[3];
sx q[3];
rz(-1.1465237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65761143) q[0];
sx q[0];
rz(-2.681499) q[0];
sx q[0];
rz(-2.7125603) q[0];
rz(1.8589164) q[1];
sx q[1];
rz(-1.7762643) q[1];
sx q[1];
rz(-0.081258953) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0690828) q[0];
sx q[0];
rz(-1.4920425) q[0];
sx q[0];
rz(2.9722636) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1090515) q[2];
sx q[2];
rz(-1.1109118) q[2];
sx q[2];
rz(2.3655917) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9489021) q[1];
sx q[1];
rz(-1.9771426) q[1];
sx q[1];
rz(-0.82273107) q[1];
rz(2.0753383) q[3];
sx q[3];
rz(-1.9849376) q[3];
sx q[3];
rz(2.4676691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.73605865) q[2];
sx q[2];
rz(-2.2106876) q[2];
sx q[2];
rz(2.0299358) q[2];
rz(-1.9723655) q[3];
sx q[3];
rz(-2.4284095) q[3];
sx q[3];
rz(-2.9986103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1375167) q[0];
sx q[0];
rz(-2.5181894) q[0];
sx q[0];
rz(-0.85859437) q[0];
rz(0.89523347) q[1];
sx q[1];
rz(-1.3139071) q[1];
sx q[1];
rz(0.039176686) q[1];
rz(1.8276385) q[2];
sx q[2];
rz(-1.5594635) q[2];
sx q[2];
rz(2.9008627) q[2];
rz(-2.9650727) q[3];
sx q[3];
rz(-0.64809496) q[3];
sx q[3];
rz(2.5213836) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
