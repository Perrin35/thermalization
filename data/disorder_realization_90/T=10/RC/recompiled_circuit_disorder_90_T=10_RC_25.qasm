OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4164299) q[0];
sx q[0];
rz(-0.13983146) q[0];
sx q[0];
rz(-2.5319985) q[0];
rz(0.66863376) q[1];
sx q[1];
rz(-2.2761087) q[1];
sx q[1];
rz(-0.087021526) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5370109) q[0];
sx q[0];
rz(-0.7365948) q[0];
sx q[0];
rz(2.1405311) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6992703) q[2];
sx q[2];
rz(-0.80768425) q[2];
sx q[2];
rz(1.8345923) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7510208) q[1];
sx q[1];
rz(-2.6387072) q[1];
sx q[1];
rz(1.2807756) q[1];
x q[2];
rz(-1.5243516) q[3];
sx q[3];
rz(-1.159045) q[3];
sx q[3];
rz(2.9573033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0027851) q[2];
sx q[2];
rz(-1.7613208) q[2];
sx q[2];
rz(0.37386093) q[2];
rz(-0.3368245) q[3];
sx q[3];
rz(-1.5954433) q[3];
sx q[3];
rz(-0.22836223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8841298) q[0];
sx q[0];
rz(-0.3733491) q[0];
sx q[0];
rz(1.194838) q[0];
rz(0.082611235) q[1];
sx q[1];
rz(-1.1673085) q[1];
sx q[1];
rz(3.1412178) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21391348) q[0];
sx q[0];
rz(-1.1535026) q[0];
sx q[0];
rz(0.37950619) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2194901) q[2];
sx q[2];
rz(-0.32425913) q[2];
sx q[2];
rz(2.1622554) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0034069) q[1];
sx q[1];
rz(-0.97461838) q[1];
sx q[1];
rz(-1.2852933) q[1];
x q[2];
rz(1.5156636) q[3];
sx q[3];
rz(-2.030636) q[3];
sx q[3];
rz(0.22051792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.80883819) q[2];
sx q[2];
rz(-0.19503441) q[2];
sx q[2];
rz(0.17671281) q[2];
rz(2.3475032) q[3];
sx q[3];
rz(-2.3637171) q[3];
sx q[3];
rz(2.5879522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2394543) q[0];
sx q[0];
rz(-0.70394009) q[0];
sx q[0];
rz(2.7368271) q[0];
rz(1.2813214) q[1];
sx q[1];
rz(-0.57360137) q[1];
sx q[1];
rz(-1.8331029) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8264016) q[0];
sx q[0];
rz(-2.0668525) q[0];
sx q[0];
rz(-1.0415003) q[0];
rz(-pi) q[1];
rz(-1.3345509) q[2];
sx q[2];
rz(-1.0457195) q[2];
sx q[2];
rz(-0.46307785) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1723459) q[1];
sx q[1];
rz(-1.8241276) q[1];
sx q[1];
rz(-2.2971056) q[1];
rz(-pi) q[2];
rz(0.19034068) q[3];
sx q[3];
rz(-2.6958709) q[3];
sx q[3];
rz(2.1325071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9222766) q[2];
sx q[2];
rz(-0.52307659) q[2];
sx q[2];
rz(-0.02040872) q[2];
rz(1.071788) q[3];
sx q[3];
rz(-2.0139549) q[3];
sx q[3];
rz(-2.4782457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9999009) q[0];
sx q[0];
rz(-0.90943709) q[0];
sx q[0];
rz(-0.91598696) q[0];
rz(0.46332106) q[1];
sx q[1];
rz(-1.0659734) q[1];
sx q[1];
rz(-2.0844918) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0650478) q[0];
sx q[0];
rz(-2.0079945) q[0];
sx q[0];
rz(1.5426427) q[0];
x q[1];
rz(0.57025036) q[2];
sx q[2];
rz(-1.0409365) q[2];
sx q[2];
rz(-2.3392764) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.108236) q[1];
sx q[1];
rz(-1.5844966) q[1];
sx q[1];
rz(-1.3407047) q[1];
x q[2];
rz(-2.9017157) q[3];
sx q[3];
rz(-2.3696218) q[3];
sx q[3];
rz(-0.039226942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.89381924) q[2];
sx q[2];
rz(-0.63988581) q[2];
sx q[2];
rz(-0.34269732) q[2];
rz(-1.4438859) q[3];
sx q[3];
rz(-0.87564898) q[3];
sx q[3];
rz(0.7152344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7856359) q[0];
sx q[0];
rz(-2.5781093) q[0];
sx q[0];
rz(0.086439565) q[0];
rz(1.3899639) q[1];
sx q[1];
rz(-0.93170634) q[1];
sx q[1];
rz(-0.46868971) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2738004) q[0];
sx q[0];
rz(-2.3228354) q[0];
sx q[0];
rz(1.0206945) q[0];
rz(-pi) q[1];
rz(-2.6195171) q[2];
sx q[2];
rz(-2.7895088) q[2];
sx q[2];
rz(0.075369518) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9710755) q[1];
sx q[1];
rz(-1.3689539) q[1];
sx q[1];
rz(-1.7032196) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0607871) q[3];
sx q[3];
rz(-0.81597933) q[3];
sx q[3];
rz(-0.24085837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9995352) q[2];
sx q[2];
rz(-2.2613328) q[2];
sx q[2];
rz(-2.2772677) q[2];
rz(2.9243829) q[3];
sx q[3];
rz(-2.5253798) q[3];
sx q[3];
rz(-0.29278452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42751673) q[0];
sx q[0];
rz(-1.7345411) q[0];
sx q[0];
rz(2.9445904) q[0];
rz(-1.3621832) q[1];
sx q[1];
rz(-1.660659) q[1];
sx q[1];
rz(0.33624712) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2339736) q[0];
sx q[0];
rz(-2.189878) q[0];
sx q[0];
rz(2.1179849) q[0];
x q[1];
rz(-2.5480812) q[2];
sx q[2];
rz(-2.0715908) q[2];
sx q[2];
rz(-1.8261432) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0974478) q[1];
sx q[1];
rz(-2.1598585) q[1];
sx q[1];
rz(1.5138022) q[1];
rz(2.9365262) q[3];
sx q[3];
rz(-0.70687095) q[3];
sx q[3];
rz(2.6170078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6062935) q[2];
sx q[2];
rz(-1.4761816) q[2];
sx q[2];
rz(-2.2952648) q[2];
rz(-1.2396631) q[3];
sx q[3];
rz(-1.9097493) q[3];
sx q[3];
rz(0.0011750778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7476615) q[0];
sx q[0];
rz(-2.1037536) q[0];
sx q[0];
rz(0.25892648) q[0];
rz(1.3461643) q[1];
sx q[1];
rz(-1.3849473) q[1];
sx q[1];
rz(-2.0475725) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0712229) q[0];
sx q[0];
rz(-0.56814146) q[0];
sx q[0];
rz(0.90026654) q[0];
rz(-0.012374087) q[2];
sx q[2];
rz(-1.3482598) q[2];
sx q[2];
rz(0.48977938) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0371857) q[1];
sx q[1];
rz(-1.6652602) q[1];
sx q[1];
rz(3.1135998) q[1];
x q[2];
rz(1.592698) q[3];
sx q[3];
rz(-1.6688445) q[3];
sx q[3];
rz(-2.736562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9672433) q[2];
sx q[2];
rz(-1.3568342) q[2];
sx q[2];
rz(2.8209177) q[2];
rz(-0.43618068) q[3];
sx q[3];
rz(-0.47302055) q[3];
sx q[3];
rz(0.44803739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86673474) q[0];
sx q[0];
rz(-0.47645706) q[0];
sx q[0];
rz(2.136769) q[0];
rz(0.59016219) q[1];
sx q[1];
rz(-0.90548038) q[1];
sx q[1];
rz(2.337713) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.852254) q[0];
sx q[0];
rz(-0.58710557) q[0];
sx q[0];
rz(3.0703074) q[0];
x q[1];
rz(-1.7558394) q[2];
sx q[2];
rz(-1.0244601) q[2];
sx q[2];
rz(-1.1748558) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4623973) q[1];
sx q[1];
rz(-3.052366) q[1];
sx q[1];
rz(1.6855082) q[1];
x q[2];
rz(1.6611093) q[3];
sx q[3];
rz(-2.2323425) q[3];
sx q[3];
rz(0.46935287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8089495) q[2];
sx q[2];
rz(-0.84471622) q[2];
sx q[2];
rz(-2.1311029) q[2];
rz(0.60339749) q[3];
sx q[3];
rz(-1.5248652) q[3];
sx q[3];
rz(2.5914014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6036966) q[0];
sx q[0];
rz(-1.864707) q[0];
sx q[0];
rz(-2.7340775) q[0];
rz(-2.852476) q[1];
sx q[1];
rz(-2.018785) q[1];
sx q[1];
rz(-2.3908652) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6347046) q[0];
sx q[0];
rz(-1.2872211) q[0];
sx q[0];
rz(2.8371235) q[0];
rz(-pi) q[1];
rz(-0.27751343) q[2];
sx q[2];
rz(-0.29007426) q[2];
sx q[2];
rz(1.8857423) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1474485) q[1];
sx q[1];
rz(-1.3683043) q[1];
sx q[1];
rz(1.2575498) q[1];
x q[2];
rz(1.5446072) q[3];
sx q[3];
rz(-1.2695754) q[3];
sx q[3];
rz(1.6645886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0687381) q[2];
sx q[2];
rz(-0.72454238) q[2];
sx q[2];
rz(-1.9753974) q[2];
rz(1.7769622) q[3];
sx q[3];
rz(-0.77562538) q[3];
sx q[3];
rz(-0.021818074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5584548) q[0];
sx q[0];
rz(-0.82644176) q[0];
sx q[0];
rz(1.3903842) q[0];
rz(-2.8109) q[1];
sx q[1];
rz(-0.76534098) q[1];
sx q[1];
rz(1.6814544) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5064142) q[0];
sx q[0];
rz(-1.9902475) q[0];
sx q[0];
rz(2.123453) q[0];
rz(-1.6129458) q[2];
sx q[2];
rz(-2.6372006) q[2];
sx q[2];
rz(-1.991589) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.32669386) q[1];
sx q[1];
rz(-1.1732475) q[1];
sx q[1];
rz(1.3709929) q[1];
rz(-pi) q[2];
rz(0.71698935) q[3];
sx q[3];
rz(-2.2319712) q[3];
sx q[3];
rz(-0.17718525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.34974393) q[2];
sx q[2];
rz(-1.8169553) q[2];
sx q[2];
rz(1.6798518) q[2];
rz(2.0215624) q[3];
sx q[3];
rz(-0.62168613) q[3];
sx q[3];
rz(-0.45564836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5158952) q[0];
sx q[0];
rz(-1.5213756) q[0];
sx q[0];
rz(-1.9949927) q[0];
rz(-1.760578) q[1];
sx q[1];
rz(-1.2881423) q[1];
sx q[1];
rz(-1.2013411) q[1];
rz(1.21576) q[2];
sx q[2];
rz(-2.3264865) q[2];
sx q[2];
rz(-2.9441499) q[2];
rz(1.1827042) q[3];
sx q[3];
rz(-2.2065065) q[3];
sx q[3];
rz(1.828215) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];