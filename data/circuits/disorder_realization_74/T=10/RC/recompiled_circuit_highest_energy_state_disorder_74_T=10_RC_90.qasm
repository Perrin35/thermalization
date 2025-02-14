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
rz(1.9407152) q[0];
sx q[0];
rz(-2.241029) q[0];
sx q[0];
rz(-2.9086034) q[0];
rz(-1.5268582) q[1];
sx q[1];
rz(4.2136636) q[1];
sx q[1];
rz(16.827488) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096123556) q[0];
sx q[0];
rz(-1.5579559) q[0];
sx q[0];
rz(1.8052285) q[0];
rz(-1.5824116) q[2];
sx q[2];
rz(-1.0093952) q[2];
sx q[2];
rz(2.7593768) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0927375) q[1];
sx q[1];
rz(-1.6462781) q[1];
sx q[1];
rz(1.9239747) q[1];
rz(0.60910881) q[3];
sx q[3];
rz(-1.2909596) q[3];
sx q[3];
rz(-2.5697054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.69922525) q[2];
sx q[2];
rz(-2.0825601) q[2];
sx q[2];
rz(3.0287058) q[2];
rz(0.26116192) q[3];
sx q[3];
rz(-1.7966725) q[3];
sx q[3];
rz(-1.2507218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3498822) q[0];
sx q[0];
rz(-0.078502027) q[0];
sx q[0];
rz(-3.0872524) q[0];
rz(2.9229274) q[1];
sx q[1];
rz(-1.668914) q[1];
sx q[1];
rz(0.36453077) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35848356) q[0];
sx q[0];
rz(-0.82573422) q[0];
sx q[0];
rz(-2.4381105) q[0];
x q[1];
rz(-2.7414315) q[2];
sx q[2];
rz(-2.1605943) q[2];
sx q[2];
rz(1.1584182) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6682525) q[1];
sx q[1];
rz(-2.1917731) q[1];
sx q[1];
rz(1.4771858) q[1];
x q[2];
rz(1.0899312) q[3];
sx q[3];
rz(-0.37853795) q[3];
sx q[3];
rz(2.7906281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.24078044) q[2];
sx q[2];
rz(-1.4419) q[2];
sx q[2];
rz(-2.0330632) q[2];
rz(2.945914) q[3];
sx q[3];
rz(-1.9340197) q[3];
sx q[3];
rz(-1.4669363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-2.7473258) q[0];
sx q[0];
rz(-1.0402004) q[0];
sx q[0];
rz(-1.0362097) q[0];
rz(-0.58468435) q[1];
sx q[1];
rz(-1.4671289) q[1];
sx q[1];
rz(2.5856957) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3506539) q[0];
sx q[0];
rz(-2.2430111) q[0];
sx q[0];
rz(1.0291128) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4934718) q[2];
sx q[2];
rz(-1.5333042) q[2];
sx q[2];
rz(0.24399569) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.039274065) q[1];
sx q[1];
rz(-2.2688451) q[1];
sx q[1];
rz(-0.1711425) q[1];
x q[2];
rz(-0.44938748) q[3];
sx q[3];
rz(-1.3419749) q[3];
sx q[3];
rz(-0.23923161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7613775) q[2];
sx q[2];
rz(-0.20647241) q[2];
sx q[2];
rz(-1.192344) q[2];
rz(-3.0377667) q[3];
sx q[3];
rz(-1.1622279) q[3];
sx q[3];
rz(1.1473568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.035606774) q[0];
sx q[0];
rz(-2.5150531) q[0];
sx q[0];
rz(-1.71738) q[0];
rz(0.72319889) q[1];
sx q[1];
rz(-0.23591787) q[1];
sx q[1];
rz(-0.22044388) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.854326) q[0];
sx q[0];
rz(-2.0279071) q[0];
sx q[0];
rz(-0.4099538) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29517572) q[2];
sx q[2];
rz(-2.1116426) q[2];
sx q[2];
rz(1.6991655) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.38762906) q[1];
sx q[1];
rz(-2.7668608) q[1];
sx q[1];
rz(1.8829569) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7679522) q[3];
sx q[3];
rz(-1.9585525) q[3];
sx q[3];
rz(1.0755634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.38349884) q[2];
sx q[2];
rz(-1.3145964) q[2];
sx q[2];
rz(-2.626075) q[2];
rz(3.0564485) q[3];
sx q[3];
rz(-2.4862423) q[3];
sx q[3];
rz(-2.3084124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4715217) q[0];
sx q[0];
rz(-0.18297289) q[0];
sx q[0];
rz(0.66437379) q[0];
rz(2.6145256) q[1];
sx q[1];
rz(-2.0030237) q[1];
sx q[1];
rz(-1.1767496) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4445317) q[0];
sx q[0];
rz(-0.089655487) q[0];
sx q[0];
rz(1.2517002) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88567733) q[2];
sx q[2];
rz(-1.9373477) q[2];
sx q[2];
rz(-1.7358071) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.47065319) q[1];
sx q[1];
rz(-1.3325508) q[1];
sx q[1];
rz(1.5907423) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99915766) q[3];
sx q[3];
rz(-1.7222341) q[3];
sx q[3];
rz(-3.0157523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.027017) q[2];
sx q[2];
rz(-2.9288374) q[2];
sx q[2];
rz(2.2034755) q[2];
rz(1.045687) q[3];
sx q[3];
rz(-2.9595879) q[3];
sx q[3];
rz(-2.3312881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.806458) q[0];
sx q[0];
rz(-1.9430176) q[0];
sx q[0];
rz(1.8916116) q[0];
rz(-2.9373923) q[1];
sx q[1];
rz(-0.67664346) q[1];
sx q[1];
rz(0.85320371) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6224504) q[0];
sx q[0];
rz(-1.2850437) q[0];
sx q[0];
rz(-2.7842872) q[0];
rz(-2.3647652) q[2];
sx q[2];
rz(-2.3558801) q[2];
sx q[2];
rz(1.7091027) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0895929) q[1];
sx q[1];
rz(-1.1530563) q[1];
sx q[1];
rz(-0.79302782) q[1];
rz(-3.0141287) q[3];
sx q[3];
rz(-1.9098305) q[3];
sx q[3];
rz(1.6515428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0975254) q[2];
sx q[2];
rz(-1.3419469) q[2];
sx q[2];
rz(1.1996783) q[2];
rz(-1.0058282) q[3];
sx q[3];
rz(-0.89322105) q[3];
sx q[3];
rz(-0.83542663) q[3];
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
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6658039) q[0];
sx q[0];
rz(-2.8732193) q[0];
sx q[0];
rz(2.1589808) q[0];
rz(-1.8334552) q[1];
sx q[1];
rz(-1.9255226) q[1];
sx q[1];
rz(-1.7024202) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7076232) q[0];
sx q[0];
rz(-0.48309193) q[0];
sx q[0];
rz(-0.77204319) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7050883) q[2];
sx q[2];
rz(-0.52596131) q[2];
sx q[2];
rz(0.44912042) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8285259) q[1];
sx q[1];
rz(-1.2432352) q[1];
sx q[1];
rz(0.86159535) q[1];
rz(-pi) q[2];
rz(-1.7129219) q[3];
sx q[3];
rz(-1.097659) q[3];
sx q[3];
rz(-2.9700235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.16284379) q[2];
sx q[2];
rz(-0.89357251) q[2];
sx q[2];
rz(2.5992375) q[2];
rz(-2.1584623) q[3];
sx q[3];
rz(-1.977908) q[3];
sx q[3];
rz(-0.94016176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5777609) q[0];
sx q[0];
rz(-1.7153808) q[0];
sx q[0];
rz(-0.78052178) q[0];
rz(1.1753987) q[1];
sx q[1];
rz(-0.54882097) q[1];
sx q[1];
rz(-1.0343879) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8261125) q[0];
sx q[0];
rz(-1.9404991) q[0];
sx q[0];
rz(-1.3913054) q[0];
x q[1];
rz(-2.0883191) q[2];
sx q[2];
rz(-0.1642326) q[2];
sx q[2];
rz(2.9688778) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.24496291) q[1];
sx q[1];
rz(-1.1707998) q[1];
sx q[1];
rz(1.4463615) q[1];
rz(-0.54776056) q[3];
sx q[3];
rz(-1.3410853) q[3];
sx q[3];
rz(-1.821777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.15022755) q[2];
sx q[2];
rz(-1.4069724) q[2];
sx q[2];
rz(0.11037174) q[2];
rz(-1.8433833) q[3];
sx q[3];
rz(-1.7896264) q[3];
sx q[3];
rz(2.8492294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89123911) q[0];
sx q[0];
rz(-1.013101) q[0];
sx q[0];
rz(1.2497586) q[0];
rz(0.091014422) q[1];
sx q[1];
rz(-2.4200771) q[1];
sx q[1];
rz(-2.4299842) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34214333) q[0];
sx q[0];
rz(-2.3490688) q[0];
sx q[0];
rz(-1.5802797) q[0];
rz(-2.4170824) q[2];
sx q[2];
rz(-2.2045928) q[2];
sx q[2];
rz(-0.29622918) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.27798692) q[1];
sx q[1];
rz(-2.5558439) q[1];
sx q[1];
rz(0.65237712) q[1];
x q[2];
rz(2.3859714) q[3];
sx q[3];
rz(-1.8241624) q[3];
sx q[3];
rz(-0.33035914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8083501) q[2];
sx q[2];
rz(-1.4834206) q[2];
sx q[2];
rz(-1.1727772) q[2];
rz(-1.6138389) q[3];
sx q[3];
rz(-2.0634191) q[3];
sx q[3];
rz(3.0465904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.461819) q[0];
sx q[0];
rz(-3.0175896) q[0];
sx q[0];
rz(1.5754196) q[0];
rz(0.35789403) q[1];
sx q[1];
rz(-1.0708258) q[1];
sx q[1];
rz(2.2391052) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76644635) q[0];
sx q[0];
rz(-1.3129108) q[0];
sx q[0];
rz(2.2058486) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4975791) q[2];
sx q[2];
rz(-1.5531544) q[2];
sx q[2];
rz(2.0669075) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8621626) q[1];
sx q[1];
rz(-1.9050652) q[1];
sx q[1];
rz(0.80252711) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0542576) q[3];
sx q[3];
rz(-1.0221586) q[3];
sx q[3];
rz(-0.61591301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1050528) q[2];
sx q[2];
rz(-2.5532711) q[2];
sx q[2];
rz(1.6416637) q[2];
rz(2.6203652) q[3];
sx q[3];
rz(-0.14385496) q[3];
sx q[3];
rz(-0.95411333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70984107) q[0];
sx q[0];
rz(-2.3713645) q[0];
sx q[0];
rz(1.4159528) q[0];
rz(3.1406828) q[1];
sx q[1];
rz(-1.6699427) q[1];
sx q[1];
rz(-1.4572399) q[1];
rz(1.4191237) q[2];
sx q[2];
rz(-2.5628452) q[2];
sx q[2];
rz(2.645523) q[2];
rz(0.72003638) q[3];
sx q[3];
rz(-0.89012082) q[3];
sx q[3];
rz(-0.16142798) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
