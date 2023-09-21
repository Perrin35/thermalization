OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8437334) q[0];
sx q[0];
rz(-0.61367404) q[0];
sx q[0];
rz(-2.4224129) q[0];
rz(1.367388) q[1];
sx q[1];
rz(-0.24582882) q[1];
sx q[1];
rz(2.153102) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2253217) q[0];
sx q[0];
rz(-2.7164408) q[0];
sx q[0];
rz(-1.4293554) q[0];
rz(1.3349873) q[2];
sx q[2];
rz(-0.60496444) q[2];
sx q[2];
rz(-1.4750823) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.72585427) q[1];
sx q[1];
rz(-1.1032747) q[1];
sx q[1];
rz(-2.095925) q[1];
rz(1.418781) q[3];
sx q[3];
rz(-2.8957643) q[3];
sx q[3];
rz(-1.0600526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6926379) q[2];
sx q[2];
rz(-1.8779034) q[2];
sx q[2];
rz(-0.56837481) q[2];
rz(1.189399) q[3];
sx q[3];
rz(-2.9247734) q[3];
sx q[3];
rz(2.9590759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4028567) q[0];
sx q[0];
rz(-2.0500654) q[0];
sx q[0];
rz(-0.36303315) q[0];
rz(-0.96827132) q[1];
sx q[1];
rz(-2.4666511) q[1];
sx q[1];
rz(-1.8889069) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7583) q[0];
sx q[0];
rz(-2.9581666) q[0];
sx q[0];
rz(0.51390506) q[0];
x q[1];
rz(0.57931309) q[2];
sx q[2];
rz(-0.99537731) q[2];
sx q[2];
rz(-2.343611) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9892271) q[1];
sx q[1];
rz(-2.3902635) q[1];
sx q[1];
rz(0.82528021) q[1];
rz(-pi) q[2];
rz(1.1870175) q[3];
sx q[3];
rz(-0.46758258) q[3];
sx q[3];
rz(3.0083857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0210555) q[2];
sx q[2];
rz(-1.3201822) q[2];
sx q[2];
rz(-0.34376124) q[2];
rz(-2.8072642) q[3];
sx q[3];
rz(-2.7407586) q[3];
sx q[3];
rz(-0.46935579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47135982) q[0];
sx q[0];
rz(-2.8911262) q[0];
sx q[0];
rz(3.1345471) q[0];
rz(2.7650611) q[1];
sx q[1];
rz(-0.9286325) q[1];
sx q[1];
rz(0.71281707) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26582742) q[0];
sx q[0];
rz(-1.7808502) q[0];
sx q[0];
rz(1.34583) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1926646) q[2];
sx q[2];
rz(-1.4107804) q[2];
sx q[2];
rz(-0.68607054) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7506999) q[1];
sx q[1];
rz(-2.0999968) q[1];
sx q[1];
rz(-0.23551029) q[1];
rz(-1.4180693) q[3];
sx q[3];
rz(-2.6543791) q[3];
sx q[3];
rz(1.0917851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.787848) q[2];
sx q[2];
rz(-1.5793844) q[2];
sx q[2];
rz(2.4776069) q[2];
rz(1.239423) q[3];
sx q[3];
rz(-0.23012161) q[3];
sx q[3];
rz(-0.91528875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5749213) q[0];
sx q[0];
rz(-0.029733812) q[0];
sx q[0];
rz(-2.5675039) q[0];
rz(0.29218778) q[1];
sx q[1];
rz(-0.27583396) q[1];
sx q[1];
rz(-2.2132197) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61361109) q[0];
sx q[0];
rz(-0.70603849) q[0];
sx q[0];
rz(1.5692977) q[0];
rz(-pi) q[1];
rz(-2.8027014) q[2];
sx q[2];
rz(-0.8182943) q[2];
sx q[2];
rz(1.3202867) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.91765431) q[1];
sx q[1];
rz(-1.8173216) q[1];
sx q[1];
rz(-1.1397821) q[1];
rz(2.2286962) q[3];
sx q[3];
rz(-1.9808589) q[3];
sx q[3];
rz(-1.0421154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0488247) q[2];
sx q[2];
rz(-0.51468819) q[2];
sx q[2];
rz(1.9862004) q[2];
rz(-2.998735) q[3];
sx q[3];
rz(-0.5345878) q[3];
sx q[3];
rz(2.3891881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083754152) q[0];
sx q[0];
rz(-2.6733119) q[0];
sx q[0];
rz(0.074247867) q[0];
rz(1.6429365) q[1];
sx q[1];
rz(-1.628412) q[1];
sx q[1];
rz(-2.1898851) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.689752) q[0];
sx q[0];
rz(-1.6432439) q[0];
sx q[0];
rz(2.3206582) q[0];
x q[1];
rz(0.81163002) q[2];
sx q[2];
rz(-1.1434165) q[2];
sx q[2];
rz(2.6805834) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.73478991) q[1];
sx q[1];
rz(-2.4996335) q[1];
sx q[1];
rz(-1.2924679) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0025025) q[3];
sx q[3];
rz(-1.2270524) q[3];
sx q[3];
rz(-0.4337726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4962861) q[2];
sx q[2];
rz(-0.22660613) q[2];
sx q[2];
rz(-1.9892233) q[2];
rz(-2.0955829) q[3];
sx q[3];
rz(-2.4029616) q[3];
sx q[3];
rz(-3.1402821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-2.1482658) q[0];
sx q[0];
rz(-2.1874805) q[0];
sx q[0];
rz(-2.6519725) q[0];
rz(-2.1173677) q[1];
sx q[1];
rz(-2.5074597) q[1];
sx q[1];
rz(1.5354059) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4682966) q[0];
sx q[0];
rz(-1.6185986) q[0];
sx q[0];
rz(3.0424546) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3697482) q[2];
sx q[2];
rz(-2.2945885) q[2];
sx q[2];
rz(0.48230241) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8294576) q[1];
sx q[1];
rz(-1.2873988) q[1];
sx q[1];
rz(-2.9764922) q[1];
x q[2];
rz(-0.225004) q[3];
sx q[3];
rz(-2.5000754) q[3];
sx q[3];
rz(-1.0596421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.17807047) q[2];
sx q[2];
rz(-2.3242951) q[2];
sx q[2];
rz(2.9658588) q[2];
rz(2.0641616) q[3];
sx q[3];
rz(-1.146799) q[3];
sx q[3];
rz(-3.1350465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0734633) q[0];
sx q[0];
rz(-0.21502762) q[0];
sx q[0];
rz(-1.8716795) q[0];
rz(0.1779671) q[1];
sx q[1];
rz(-2.3033419) q[1];
sx q[1];
rz(-1.3659182) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89643909) q[0];
sx q[0];
rz(-1.23587) q[0];
sx q[0];
rz(2.4048637) q[0];
rz(-pi) q[1];
rz(0.62909158) q[2];
sx q[2];
rz(-2.7436939) q[2];
sx q[2];
rz(-1.2049904) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1508031) q[1];
sx q[1];
rz(-1.8263131) q[1];
sx q[1];
rz(-1.7533416) q[1];
rz(-2.9317022) q[3];
sx q[3];
rz(-0.53901796) q[3];
sx q[3];
rz(-0.15912661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0718677) q[2];
sx q[2];
rz(-0.29399997) q[2];
sx q[2];
rz(-2.243637) q[2];
rz(-0.14363025) q[3];
sx q[3];
rz(-1.2575905) q[3];
sx q[3];
rz(-0.34415054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1865858) q[0];
sx q[0];
rz(-0.85383677) q[0];
sx q[0];
rz(-0.32178497) q[0];
rz(2.2161662) q[1];
sx q[1];
rz(-1.7089475) q[1];
sx q[1];
rz(2.5245573) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5962284) q[0];
sx q[0];
rz(-1.3188688) q[0];
sx q[0];
rz(2.9357301) q[0];
rz(1.4204942) q[2];
sx q[2];
rz(-2.9234924) q[2];
sx q[2];
rz(0.66767603) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2136692) q[1];
sx q[1];
rz(-2.546715) q[1];
sx q[1];
rz(-2.2713714) q[1];
x q[2];
rz(1.5259605) q[3];
sx q[3];
rz(-0.82190824) q[3];
sx q[3];
rz(0.84079784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.59163219) q[2];
sx q[2];
rz(-2.154921) q[2];
sx q[2];
rz(0.11403306) q[2];
rz(0.36241254) q[3];
sx q[3];
rz(-2.896538) q[3];
sx q[3];
rz(1.2362278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47700259) q[0];
sx q[0];
rz(-0.53628558) q[0];
sx q[0];
rz(-1.7846918) q[0];
rz(0.82018954) q[1];
sx q[1];
rz(-0.35239041) q[1];
sx q[1];
rz(-1.6581416) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9634919) q[0];
sx q[0];
rz(-0.083847001) q[0];
sx q[0];
rz(1.447178) q[0];
rz(2.9786417) q[2];
sx q[2];
rz(-0.59646791) q[2];
sx q[2];
rz(2.5359254) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.22944731) q[1];
sx q[1];
rz(-1.427269) q[1];
sx q[1];
rz(1.9043282) q[1];
x q[2];
rz(-1.9547192) q[3];
sx q[3];
rz(-2.0878289) q[3];
sx q[3];
rz(0.25728713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1411529) q[2];
sx q[2];
rz(-2.4856462) q[2];
sx q[2];
rz(-1.8971987) q[2];
rz(2.7164298) q[3];
sx q[3];
rz(-2.3624698) q[3];
sx q[3];
rz(1.6311133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(2.4193831) q[0];
sx q[0];
rz(-2.6477551) q[0];
sx q[0];
rz(2.7710932) q[0];
rz(-1.1100948) q[1];
sx q[1];
rz(-0.67567527) q[1];
sx q[1];
rz(1.3409021) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6275218) q[0];
sx q[0];
rz(-0.43826575) q[0];
sx q[0];
rz(-1.0036432) q[0];
rz(-2.9843763) q[2];
sx q[2];
rz(-0.72050691) q[2];
sx q[2];
rz(0.23888982) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.51532981) q[1];
sx q[1];
rz(-1.9326107) q[1];
sx q[1];
rz(-2.8742909) q[1];
rz(-2.9756536) q[3];
sx q[3];
rz(-0.12291848) q[3];
sx q[3];
rz(-2.3636706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5122539) q[2];
sx q[2];
rz(-2.9014355) q[2];
sx q[2];
rz(-1.6188999) q[2];
rz(-2.2807138) q[3];
sx q[3];
rz(-1.7247793) q[3];
sx q[3];
rz(-0.035877429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9027949) q[0];
sx q[0];
rz(-1.4807985) q[0];
sx q[0];
rz(2.280622) q[0];
rz(-2.8181656) q[1];
sx q[1];
rz(-1.885781) q[1];
sx q[1];
rz(1.5992004) q[1];
rz(1.515083) q[2];
sx q[2];
rz(-0.40691661) q[2];
sx q[2];
rz(-0.42452068) q[2];
rz(-0.38701804) q[3];
sx q[3];
rz(-1.1292463) q[3];
sx q[3];
rz(-2.675727) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
