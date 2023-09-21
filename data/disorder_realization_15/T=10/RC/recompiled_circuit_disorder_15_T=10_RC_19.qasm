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
rz(0.71917978) q[0];
rz(1.367388) q[1];
sx q[1];
rz(-0.24582882) q[1];
sx q[1];
rz(2.153102) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6160625) q[0];
sx q[0];
rz(-1.6289734) q[0];
sx q[0];
rz(1.9921897) q[0];
rz(-pi) q[1];
rz(-2.9814331) q[2];
sx q[2];
rz(-2.1567492) q[2];
sx q[2];
rz(-1.9507267) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6355675) q[1];
sx q[1];
rz(-0.68817455) q[1];
sx q[1];
rz(-0.78189214) q[1];
x q[2];
rz(3.1036166) q[3];
sx q[3];
rz(-1.8137323) q[3];
sx q[3];
rz(-1.9248885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6926379) q[2];
sx q[2];
rz(-1.2636893) q[2];
sx q[2];
rz(2.5732178) q[2];
rz(-1.189399) q[3];
sx q[3];
rz(-2.9247734) q[3];
sx q[3];
rz(0.18251671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73873591) q[0];
sx q[0];
rz(-2.0500654) q[0];
sx q[0];
rz(-0.36303315) q[0];
rz(2.1733213) q[1];
sx q[1];
rz(-0.67494154) q[1];
sx q[1];
rz(-1.2526858) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68080901) q[0];
sx q[0];
rz(-1.4810116) q[0];
sx q[0];
rz(0.16016527) q[0];
rz(0.86979903) q[2];
sx q[2];
rz(-2.3491078) q[2];
sx q[2];
rz(-3.062641) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.82488793) q[1];
sx q[1];
rz(-2.0522293) q[1];
sx q[1];
rz(2.1722721) q[1];
rz(-pi) q[2];
rz(0.18685347) q[3];
sx q[3];
rz(-1.1396176) q[3];
sx q[3];
rz(-2.8499875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0210555) q[2];
sx q[2];
rz(-1.8214104) q[2];
sx q[2];
rz(-0.34376124) q[2];
rz(0.3343285) q[3];
sx q[3];
rz(-2.7407586) q[3];
sx q[3];
rz(2.6722369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
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
rz(2.6702328) q[0];
sx q[0];
rz(-2.8911262) q[0];
sx q[0];
rz(0.0070455889) q[0];
rz(-0.37653157) q[1];
sx q[1];
rz(-2.2129602) q[1];
sx q[1];
rz(-0.71281707) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56601554) q[0];
sx q[0];
rz(-0.30656719) q[0];
sx q[0];
rz(-0.80802877) q[0];
rz(-pi) q[1];
rz(1.1926646) q[2];
sx q[2];
rz(-1.4107804) q[2];
sx q[2];
rz(-2.4555221) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7506999) q[1];
sx q[1];
rz(-1.0415959) q[1];
sx q[1];
rz(2.9060824) q[1];
rz(-1.4180693) q[3];
sx q[3];
rz(-2.6543791) q[3];
sx q[3];
rz(-2.0498076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.787848) q[2];
sx q[2];
rz(-1.5622082) q[2];
sx q[2];
rz(-0.66398579) q[2];
rz(1.239423) q[3];
sx q[3];
rz(-2.911471) q[3];
sx q[3];
rz(0.91528875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5666714) q[0];
sx q[0];
rz(-0.029733812) q[0];
sx q[0];
rz(2.5675039) q[0];
rz(0.29218778) q[1];
sx q[1];
rz(-2.8657587) q[1];
sx q[1];
rz(-0.92837292) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5279816) q[0];
sx q[0];
rz(-0.70603849) q[0];
sx q[0];
rz(-1.5722949) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9119772) q[2];
sx q[2];
rz(-2.3301635) q[2];
sx q[2];
rz(2.2974643) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.54138994) q[1];
sx q[1];
rz(-1.9879513) q[1];
sx q[1];
rz(0.27020176) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9527693) q[3];
sx q[3];
rz(-2.3828155) q[3];
sx q[3];
rz(0.052386802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.092768) q[2];
sx q[2];
rz(-0.51468819) q[2];
sx q[2];
rz(-1.1553923) q[2];
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
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0578385) q[0];
sx q[0];
rz(-2.6733119) q[0];
sx q[0];
rz(-0.074247867) q[0];
rz(1.4986562) q[1];
sx q[1];
rz(-1.628412) q[1];
sx q[1];
rz(-0.9517076) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.041391011) q[0];
sx q[0];
rz(-0.75267422) q[0];
sx q[0];
rz(1.4647096) q[0];
rz(-pi) q[1];
rz(0.56065083) q[2];
sx q[2];
rz(-2.2477305) q[2];
sx q[2];
rz(-1.4844984) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0774827) q[1];
sx q[1];
rz(-2.1842842) q[1];
sx q[1];
rz(-2.9390099) q[1];
rz(-pi) q[2];
rz(2.1577155) q[3];
sx q[3];
rz(-0.65423274) q[3];
sx q[3];
rz(-0.65140843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4962861) q[2];
sx q[2];
rz(-0.22660613) q[2];
sx q[2];
rz(-1.1523694) q[2];
rz(2.0955829) q[3];
sx q[3];
rz(-2.4029616) q[3];
sx q[3];
rz(-0.0013105198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9933269) q[0];
sx q[0];
rz(-2.1874805) q[0];
sx q[0];
rz(-0.48962012) q[0];
rz(-2.1173677) q[1];
sx q[1];
rz(-0.63413292) q[1];
sx q[1];
rz(-1.5354059) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0438457) q[0];
sx q[0];
rz(-1.6698208) q[0];
sx q[0];
rz(1.6188341) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3697482) q[2];
sx q[2];
rz(-0.8470042) q[2];
sx q[2];
rz(0.48230241) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9163497) q[1];
sx q[1];
rz(-0.32685977) q[1];
sx q[1];
rz(-1.0570231) q[1];
rz(-0.225004) q[3];
sx q[3];
rz(-2.5000754) q[3];
sx q[3];
rz(2.0819506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.17807047) q[2];
sx q[2];
rz(-2.3242951) q[2];
sx q[2];
rz(-2.9658588) q[2];
rz(2.0641616) q[3];
sx q[3];
rz(-1.146799) q[3];
sx q[3];
rz(-3.1350465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068129383) q[0];
sx q[0];
rz(-0.21502762) q[0];
sx q[0];
rz(-1.8716795) q[0];
rz(0.1779671) q[1];
sx q[1];
rz(-0.83825076) q[1];
sx q[1];
rz(-1.7756745) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89643909) q[0];
sx q[0];
rz(-1.9057227) q[0];
sx q[0];
rz(2.4048637) q[0];
rz(-pi) q[1];
rz(1.3283417) q[2];
sx q[2];
rz(-1.2521267) q[2];
sx q[2];
rz(-2.6048425) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7682225) q[1];
sx q[1];
rz(-1.3942413) q[1];
sx q[1];
rz(-0.25964398) q[1];
rz(-2.6122983) q[3];
sx q[3];
rz(-1.4636453) q[3];
sx q[3];
rz(1.2308434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0718677) q[2];
sx q[2];
rz(-0.29399997) q[2];
sx q[2];
rz(-0.89795566) q[2];
rz(-0.14363025) q[3];
sx q[3];
rz(-1.2575905) q[3];
sx q[3];
rz(-0.34415054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1865858) q[0];
sx q[0];
rz(-0.85383677) q[0];
sx q[0];
rz(0.32178497) q[0];
rz(-0.92542648) q[1];
sx q[1];
rz(-1.4326452) q[1];
sx q[1];
rz(-2.5245573) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.077438146) q[0];
sx q[0];
rz(-1.371521) q[0];
sx q[0];
rz(1.827924) q[0];
x q[1];
rz(1.7865137) q[2];
sx q[2];
rz(-1.6032013) q[2];
sx q[2];
rz(-2.0916794) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1083793) q[1];
sx q[1];
rz(-1.2011659) q[1];
sx q[1];
rz(-2.0481678) q[1];
rz(-0.7493895) q[3];
sx q[3];
rz(-1.6036311) q[3];
sx q[3];
rz(2.4421305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.59163219) q[2];
sx q[2];
rz(-0.98667163) q[2];
sx q[2];
rz(-0.11403306) q[2];
rz(-2.7791801) q[3];
sx q[3];
rz(-0.24505469) q[3];
sx q[3];
rz(-1.2362278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47700259) q[0];
sx q[0];
rz(-0.53628558) q[0];
sx q[0];
rz(1.3569008) q[0];
rz(2.3214031) q[1];
sx q[1];
rz(-0.35239041) q[1];
sx q[1];
rz(-1.483451) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8394422) q[0];
sx q[0];
rz(-1.654002) q[0];
sx q[0];
rz(-3.1312301) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16295095) q[2];
sx q[2];
rz(-0.59646791) q[2];
sx q[2];
rz(0.60566723) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8497613) q[1];
sx q[1];
rz(-1.9007678) q[1];
sx q[1];
rz(2.989819) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55012196) q[3];
sx q[3];
rz(-1.9024444) q[3];
sx q[3];
rz(1.5105997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1411529) q[2];
sx q[2];
rz(-0.65594643) q[2];
sx q[2];
rz(-1.2443939) q[2];
rz(0.42516285) q[3];
sx q[3];
rz(-2.3624698) q[3];
sx q[3];
rz(-1.6311133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4193831) q[0];
sx q[0];
rz(-2.6477551) q[0];
sx q[0];
rz(0.37049946) q[0];
rz(1.1100948) q[1];
sx q[1];
rz(-0.67567527) q[1];
sx q[1];
rz(-1.3409021) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6751624) q[0];
sx q[0];
rz(-1.8008045) q[0];
sx q[0];
rz(1.1943597) q[0];
rz(0.15721639) q[2];
sx q[2];
rz(-0.72050691) q[2];
sx q[2];
rz(-2.9027028) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.51532981) q[1];
sx q[1];
rz(-1.9326107) q[1];
sx q[1];
rz(-0.26730178) q[1];
rz(2.9756536) q[3];
sx q[3];
rz(-3.0186742) q[3];
sx q[3];
rz(0.77792203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5122539) q[2];
sx q[2];
rz(-2.9014355) q[2];
sx q[2];
rz(1.5226927) q[2];
rz(0.86087888) q[3];
sx q[3];
rz(-1.4168134) q[3];
sx q[3];
rz(-3.1057152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9027949) q[0];
sx q[0];
rz(-1.6607941) q[0];
sx q[0];
rz(-0.86097062) q[0];
rz(0.32342708) q[1];
sx q[1];
rz(-1.885781) q[1];
sx q[1];
rz(1.5992004) q[1];
rz(1.6265097) q[2];
sx q[2];
rz(-2.734676) q[2];
sx q[2];
rz(2.717072) q[2];
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