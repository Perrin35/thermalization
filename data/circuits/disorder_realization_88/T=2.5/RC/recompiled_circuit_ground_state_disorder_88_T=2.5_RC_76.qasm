OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2849046) q[0];
sx q[0];
rz(-1.8104799) q[0];
sx q[0];
rz(2.3321505) q[0];
rz(-2.5419905) q[1];
sx q[1];
rz(2.0166346) q[1];
sx q[1];
rz(12.044608) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7358462) q[0];
sx q[0];
rz(-0.1400811) q[0];
sx q[0];
rz(-1.5625885) q[0];
rz(-pi) q[1];
rz(-0.56057616) q[2];
sx q[2];
rz(-1.3180705) q[2];
sx q[2];
rz(-1.2863024) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0787766) q[1];
sx q[1];
rz(-2.605781) q[1];
sx q[1];
rz(2.1242555) q[1];
rz(-2.605569) q[3];
sx q[3];
rz(-1.6451577) q[3];
sx q[3];
rz(-1.0413359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1521505) q[2];
sx q[2];
rz(-0.70789727) q[2];
sx q[2];
rz(-1.36261) q[2];
rz(-1.4710434) q[3];
sx q[3];
rz(-1.5971284) q[3];
sx q[3];
rz(0.41489261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3926369) q[0];
sx q[0];
rz(-1.3480659) q[0];
sx q[0];
rz(-2.0718527) q[0];
rz(-2.3906129) q[1];
sx q[1];
rz(-2.2396478) q[1];
sx q[1];
rz(-2.353207) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5699002) q[0];
sx q[0];
rz(-2.4909752) q[0];
sx q[0];
rz(-2.8293737) q[0];
x q[1];
rz(-0.66547243) q[2];
sx q[2];
rz(-1.62684) q[2];
sx q[2];
rz(-1.9201345) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.091011) q[1];
sx q[1];
rz(-1.7643223) q[1];
sx q[1];
rz(1.0440318) q[1];
rz(-pi) q[2];
rz(-0.69425868) q[3];
sx q[3];
rz(-0.63000796) q[3];
sx q[3];
rz(0.97351551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1145733) q[2];
sx q[2];
rz(-0.35832778) q[2];
sx q[2];
rz(2.2349854) q[2];
rz(1.00057) q[3];
sx q[3];
rz(-1.118719) q[3];
sx q[3];
rz(-0.51893836) q[3];
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
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24666102) q[0];
sx q[0];
rz(-0.38235679) q[0];
sx q[0];
rz(-1.83778) q[0];
rz(1.9815824) q[1];
sx q[1];
rz(-1.8142533) q[1];
sx q[1];
rz(-1.4132285) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0864058) q[0];
sx q[0];
rz(-0.16567437) q[0];
sx q[0];
rz(2.590986) q[0];
rz(-0.5858461) q[2];
sx q[2];
rz(-2.4163462) q[2];
sx q[2];
rz(2.9278529) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4354463) q[1];
sx q[1];
rz(-2.6455952) q[1];
sx q[1];
rz(-2.7643327) q[1];
x q[2];
rz(0.12574844) q[3];
sx q[3];
rz(-0.85966051) q[3];
sx q[3];
rz(-1.2918732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.383519) q[2];
sx q[2];
rz(-1.5531837) q[2];
sx q[2];
rz(-0.12913945) q[2];
rz(2.6376851) q[3];
sx q[3];
rz(-1.4174856) q[3];
sx q[3];
rz(-2.5655139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1103766) q[0];
sx q[0];
rz(-1.2126558) q[0];
sx q[0];
rz(0.84621286) q[0];
rz(-0.57202488) q[1];
sx q[1];
rz(-1.6366448) q[1];
sx q[1];
rz(-1.9073073) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18674984) q[0];
sx q[0];
rz(-2.7949547) q[0];
sx q[0];
rz(0.28916547) q[0];
rz(-pi) q[1];
rz(-2.8568966) q[2];
sx q[2];
rz(-1.1276334) q[2];
sx q[2];
rz(-0.62847947) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.14951359) q[1];
sx q[1];
rz(-1.7798335) q[1];
sx q[1];
rz(-1.4884218) q[1];
rz(2.3448496) q[3];
sx q[3];
rz(-1.1827743) q[3];
sx q[3];
rz(-1.4772367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7531551) q[2];
sx q[2];
rz(-2.2074102) q[2];
sx q[2];
rz(-1.4446806) q[2];
rz(1.2706903) q[3];
sx q[3];
rz(-1.5710187) q[3];
sx q[3];
rz(1.1893893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63950771) q[0];
sx q[0];
rz(-1.4690228) q[0];
sx q[0];
rz(1.0953267) q[0];
rz(0.0630088) q[1];
sx q[1];
rz(-1.6687702) q[1];
sx q[1];
rz(-2.1953886) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.584884) q[0];
sx q[0];
rz(-2.1153304) q[0];
sx q[0];
rz(3.1088215) q[0];
rz(-pi) q[1];
rz(-0.63615648) q[2];
sx q[2];
rz(-1.8811748) q[2];
sx q[2];
rz(0.39575044) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7530147) q[1];
sx q[1];
rz(-1.5756956) q[1];
sx q[1];
rz(-2.6304662) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9619476) q[3];
sx q[3];
rz(-1.7877276) q[3];
sx q[3];
rz(2.2892078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.63518628) q[2];
sx q[2];
rz(-1.0794285) q[2];
sx q[2];
rz(1.6335454) q[2];
rz(-0.15277282) q[3];
sx q[3];
rz(-1.2341876) q[3];
sx q[3];
rz(-2.2085786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2372811) q[0];
sx q[0];
rz(-2.5560684) q[0];
sx q[0];
rz(-3.0329419) q[0];
rz(0.12734224) q[1];
sx q[1];
rz(-1.2510977) q[1];
sx q[1];
rz(-2.2551575) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3868635) q[0];
sx q[0];
rz(-0.51966705) q[0];
sx q[0];
rz(0.033367446) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2904002) q[2];
sx q[2];
rz(-0.2125769) q[2];
sx q[2];
rz(-0.55586284) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.163609) q[1];
sx q[1];
rz(-0.22630461) q[1];
sx q[1];
rz(-1.3127021) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9720794) q[3];
sx q[3];
rz(-1.4614786) q[3];
sx q[3];
rz(-1.5284571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.844937) q[2];
sx q[2];
rz(-0.79534328) q[2];
sx q[2];
rz(-2.8080158) q[2];
rz(2.2070456) q[3];
sx q[3];
rz(-2.3881674) q[3];
sx q[3];
rz(2.2540588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6880671) q[0];
sx q[0];
rz(-1.6418566) q[0];
sx q[0];
rz(2.8330579) q[0];
rz(-2.535179) q[1];
sx q[1];
rz(-2.2589222) q[1];
sx q[1];
rz(-0.44953129) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7768984) q[0];
sx q[0];
rz(-2.6230271) q[0];
sx q[0];
rz(-3.0125951) q[0];
rz(-2.0386731) q[2];
sx q[2];
rz(-0.84005594) q[2];
sx q[2];
rz(0.30147027) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7923741) q[1];
sx q[1];
rz(-1.6728396) q[1];
sx q[1];
rz(0.95178638) q[1];
rz(0.94664871) q[3];
sx q[3];
rz(-0.80393857) q[3];
sx q[3];
rz(0.80292668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.061607925) q[2];
sx q[2];
rz(-1.4915165) q[2];
sx q[2];
rz(-1.9645346) q[2];
rz(-0.70147771) q[3];
sx q[3];
rz(-0.25291118) q[3];
sx q[3];
rz(0.85566068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6845282) q[0];
sx q[0];
rz(-2.6236911) q[0];
sx q[0];
rz(1.6494226) q[0];
rz(1.0555142) q[1];
sx q[1];
rz(-1.5558473) q[1];
sx q[1];
rz(0.42748705) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2023017) q[0];
sx q[0];
rz(-0.84502506) q[0];
sx q[0];
rz(-2.2875252) q[0];
rz(-pi) q[1];
rz(-0.43085499) q[2];
sx q[2];
rz(-1.7175894) q[2];
sx q[2];
rz(2.1374747) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.44459773) q[1];
sx q[1];
rz(-2.8279404) q[1];
sx q[1];
rz(-3.0693377) q[1];
rz(-pi) q[2];
rz(-2.8264846) q[3];
sx q[3];
rz(-1.3118532) q[3];
sx q[3];
rz(-1.068598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.18275729) q[2];
sx q[2];
rz(-1.8104825) q[2];
sx q[2];
rz(-1.9367564) q[2];
rz(-0.13293535) q[3];
sx q[3];
rz(-3.0429621) q[3];
sx q[3];
rz(1.0814063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8590915) q[0];
sx q[0];
rz(-1.7318672) q[0];
sx q[0];
rz(2.6891563) q[0];
rz(-0.85894194) q[1];
sx q[1];
rz(-2.5622538) q[1];
sx q[1];
rz(1.6890242) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4005139) q[0];
sx q[0];
rz(-0.79325926) q[0];
sx q[0];
rz(-0.69852065) q[0];
rz(2.1253517) q[2];
sx q[2];
rz(-1.8565389) q[2];
sx q[2];
rz(0.72140297) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7048129) q[1];
sx q[1];
rz(-1.5680934) q[1];
sx q[1];
rz(1.0648492) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3708569) q[3];
sx q[3];
rz(-2.6477034) q[3];
sx q[3];
rz(-1.0115185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9026044) q[2];
sx q[2];
rz(-1.8391823) q[2];
sx q[2];
rz(-2.7461309) q[2];
rz(2.0579193) q[3];
sx q[3];
rz(-2.5189221) q[3];
sx q[3];
rz(1.4130939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.285242) q[0];
sx q[0];
rz(-3.0265891) q[0];
sx q[0];
rz(-1.850199) q[0];
rz(-2.3639288) q[1];
sx q[1];
rz(-1.5823369) q[1];
sx q[1];
rz(1.8596328) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6586475) q[0];
sx q[0];
rz(-1.7883669) q[0];
sx q[0];
rz(1.8013784) q[0];
rz(-pi) q[1];
rz(1.3488814) q[2];
sx q[2];
rz(-0.93524987) q[2];
sx q[2];
rz(-1.9161621) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8429148) q[1];
sx q[1];
rz(-0.58393541) q[1];
sx q[1];
rz(-2.7177277) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3814209) q[3];
sx q[3];
rz(-2.6145589) q[3];
sx q[3];
rz(2.1552212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.97734863) q[2];
sx q[2];
rz(-2.7898596) q[2];
sx q[2];
rz(0.39883167) q[2];
rz(2.2073958) q[3];
sx q[3];
rz(-2.4103006) q[3];
sx q[3];
rz(-0.092223316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4451404) q[0];
sx q[0];
rz(-0.91309375) q[0];
sx q[0];
rz(0.0050807411) q[0];
rz(0.65771163) q[1];
sx q[1];
rz(-1.4124159) q[1];
sx q[1];
rz(1.1884069) q[1];
rz(-0.80398139) q[2];
sx q[2];
rz(-1.3406499) q[2];
sx q[2];
rz(-1.7743877) q[2];
rz(-0.98500959) q[3];
sx q[3];
rz(-1.5187217) q[3];
sx q[3];
rz(0.050654618) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
