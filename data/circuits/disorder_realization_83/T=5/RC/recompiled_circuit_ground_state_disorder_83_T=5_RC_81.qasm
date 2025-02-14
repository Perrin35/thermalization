OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2228425) q[0];
sx q[0];
rz(-1.9434682) q[0];
sx q[0];
rz(-0.54984468) q[0];
rz(3.0749908) q[1];
sx q[1];
rz(-1.9220592) q[1];
sx q[1];
rz(-1.2256844) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9860814) q[0];
sx q[0];
rz(-1.5973463) q[0];
sx q[0];
rz(2.9355713) q[0];
rz(-pi) q[1];
rz(-0.066402175) q[2];
sx q[2];
rz(-1.4208671) q[2];
sx q[2];
rz(2.2663946) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5962564) q[1];
sx q[1];
rz(-1.8972555) q[1];
sx q[1];
rz(-3.1395034) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3382149) q[3];
sx q[3];
rz(-2.3720494) q[3];
sx q[3];
rz(-0.93955296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.371202) q[2];
sx q[2];
rz(-1.0545701) q[2];
sx q[2];
rz(2.2626109) q[2];
rz(0.061035872) q[3];
sx q[3];
rz(-1.6998484) q[3];
sx q[3];
rz(-0.92476168) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9691129) q[0];
sx q[0];
rz(-2.0308487) q[0];
sx q[0];
rz(-1.2193532) q[0];
rz(2.1439233) q[1];
sx q[1];
rz(-1.4439293) q[1];
sx q[1];
rz(0.77233058) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0270841) q[0];
sx q[0];
rz(-2.3128384) q[0];
sx q[0];
rz(2.7161612) q[0];
x q[1];
rz(1.9034687) q[2];
sx q[2];
rz(-1.5539697) q[2];
sx q[2];
rz(-2.1706659) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0454084) q[1];
sx q[1];
rz(-1.6077157) q[1];
sx q[1];
rz(0.30942076) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1082013) q[3];
sx q[3];
rz(-0.62761939) q[3];
sx q[3];
rz(-2.1737568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4274009) q[2];
sx q[2];
rz(-1.4560459) q[2];
sx q[2];
rz(1.6198772) q[2];
rz(1.6649668) q[3];
sx q[3];
rz(-1.0790389) q[3];
sx q[3];
rz(-2.9286706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6427226) q[0];
sx q[0];
rz(-1.9864137) q[0];
sx q[0];
rz(0.89514071) q[0];
rz(-2.8816032) q[1];
sx q[1];
rz(-1.3464059) q[1];
sx q[1];
rz(0.77494979) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4013028) q[0];
sx q[0];
rz(-2.4537114) q[0];
sx q[0];
rz(0.58346652) q[0];
rz(2.9812212) q[2];
sx q[2];
rz(-0.43284349) q[2];
sx q[2];
rz(-2.3059887) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5436478) q[1];
sx q[1];
rz(-2.0646808) q[1];
sx q[1];
rz(0.96205254) q[1];
rz(-pi) q[2];
x q[2];
rz(2.272153) q[3];
sx q[3];
rz(-2.0819391) q[3];
sx q[3];
rz(1.7366586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.2669175) q[2];
sx q[2];
rz(-0.83628925) q[2];
sx q[2];
rz(-0.8379035) q[2];
rz(-0.72495929) q[3];
sx q[3];
rz(-0.91491142) q[3];
sx q[3];
rz(-0.27423492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2570141) q[0];
sx q[0];
rz(-3.0138636) q[0];
sx q[0];
rz(1.1626441) q[0];
rz(-1.1391901) q[1];
sx q[1];
rz(-1.2290686) q[1];
sx q[1];
rz(-2.7584279) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1568147) q[0];
sx q[0];
rz(-1.1582644) q[0];
sx q[0];
rz(1.5810313) q[0];
rz(-pi) q[1];
rz(1.2960984) q[2];
sx q[2];
rz(-2.3847859) q[2];
sx q[2];
rz(-1.6414549) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.207473) q[1];
sx q[1];
rz(-1.2223244) q[1];
sx q[1];
rz(1.5847809) q[1];
rz(-2.7015436) q[3];
sx q[3];
rz(-0.69877934) q[3];
sx q[3];
rz(0.13126365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42134103) q[2];
sx q[2];
rz(-0.64695078) q[2];
sx q[2];
rz(1.7163537) q[2];
rz(1.1228115) q[3];
sx q[3];
rz(-0.69901005) q[3];
sx q[3];
rz(2.5785246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4238033) q[0];
sx q[0];
rz(-2.932982) q[0];
sx q[0];
rz(0.31722379) q[0];
rz(-0.18103655) q[1];
sx q[1];
rz(-1.92417) q[1];
sx q[1];
rz(-2.5203629) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4201403) q[0];
sx q[0];
rz(-0.83774191) q[0];
sx q[0];
rz(0.28845235) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6124837) q[2];
sx q[2];
rz(-2.0663849) q[2];
sx q[2];
rz(1.017073) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.40338755) q[1];
sx q[1];
rz(-0.59425577) q[1];
sx q[1];
rz(-1.4273093) q[1];
x q[2];
rz(-2.0192573) q[3];
sx q[3];
rz(-0.82379195) q[3];
sx q[3];
rz(2.6246678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8165596) q[2];
sx q[2];
rz(-0.72924048) q[2];
sx q[2];
rz(-1.0664335) q[2];
rz(1.019143) q[3];
sx q[3];
rz(-0.72503966) q[3];
sx q[3];
rz(-1.9371921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2524734) q[0];
sx q[0];
rz(-2.7121565) q[0];
sx q[0];
rz(-0.47527894) q[0];
rz(2.1976082) q[1];
sx q[1];
rz(-1.9119268) q[1];
sx q[1];
rz(0.42517391) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89579534) q[0];
sx q[0];
rz(-2.6962792) q[0];
sx q[0];
rz(-2.6209339) q[0];
rz(0.026907909) q[2];
sx q[2];
rz(-0.49420658) q[2];
sx q[2];
rz(1.869452) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.43229285) q[1];
sx q[1];
rz(-1.4008351) q[1];
sx q[1];
rz(-1.0372838) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75192368) q[3];
sx q[3];
rz(-2.6787191) q[3];
sx q[3];
rz(-1.183267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4117333) q[2];
sx q[2];
rz(-1.4889762) q[2];
sx q[2];
rz(0.64424166) q[2];
rz(1.1614557) q[3];
sx q[3];
rz(-2.1803653) q[3];
sx q[3];
rz(-0.78288356) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2162868) q[0];
sx q[0];
rz(-1.1290978) q[0];
sx q[0];
rz(2.570545) q[0];
rz(2.0371927) q[1];
sx q[1];
rz(-2.0490502) q[1];
sx q[1];
rz(2.8071094) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5307823) q[0];
sx q[0];
rz(-0.40662262) q[0];
sx q[0];
rz(-1.6645891) q[0];
rz(-0.38774298) q[2];
sx q[2];
rz(-0.84008145) q[2];
sx q[2];
rz(2.4304488) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.063364048) q[1];
sx q[1];
rz(-1.9029088) q[1];
sx q[1];
rz(-2.8719877) q[1];
rz(-pi) q[2];
rz(-2.7860113) q[3];
sx q[3];
rz(-0.68590859) q[3];
sx q[3];
rz(2.4816328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0337446) q[2];
sx q[2];
rz(-2.0287543) q[2];
sx q[2];
rz(-3.0372078) q[2];
rz(-2.8153822) q[3];
sx q[3];
rz(-2.2741337) q[3];
sx q[3];
rz(0.67702684) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0831083) q[0];
sx q[0];
rz(-1.9309738) q[0];
sx q[0];
rz(-2.6404358) q[0];
rz(2.0106563) q[1];
sx q[1];
rz(-2.1986745) q[1];
sx q[1];
rz(-1.2859734) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0663529) q[0];
sx q[0];
rz(-1.4202002) q[0];
sx q[0];
rz(-2.0429918) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2663009) q[2];
sx q[2];
rz(-2.2448953) q[2];
sx q[2];
rz(2.4456191) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0803538) q[1];
sx q[1];
rz(-0.29268943) q[1];
sx q[1];
rz(-0.55673843) q[1];
rz(-pi) q[2];
rz(-2.6874316) q[3];
sx q[3];
rz(-1.1060904) q[3];
sx q[3];
rz(2.2247119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8274294) q[2];
sx q[2];
rz(-0.25671998) q[2];
sx q[2];
rz(0.87744212) q[2];
rz(2.1432803) q[3];
sx q[3];
rz(-2.5613997) q[3];
sx q[3];
rz(-1.72054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5494635) q[0];
sx q[0];
rz(-1.2315467) q[0];
sx q[0];
rz(0.25729427) q[0];
rz(0.65722242) q[1];
sx q[1];
rz(-0.46923271) q[1];
sx q[1];
rz(-1.2425655) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7433194) q[0];
sx q[0];
rz(-0.87617517) q[0];
sx q[0];
rz(0.99382235) q[0];
rz(-2.9986249) q[2];
sx q[2];
rz(-0.7960628) q[2];
sx q[2];
rz(-0.23676714) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5754274) q[1];
sx q[1];
rz(-2.0708443) q[1];
sx q[1];
rz(2.0636255) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6493862) q[3];
sx q[3];
rz(-2.7653381) q[3];
sx q[3];
rz(0.87615651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2520752) q[2];
sx q[2];
rz(-1.4318848) q[2];
sx q[2];
rz(0.70402181) q[2];
rz(1.4975123) q[3];
sx q[3];
rz(-1.6234532) q[3];
sx q[3];
rz(-1.9254855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(1.7261451) q[0];
sx q[0];
rz(-0.69196597) q[0];
sx q[0];
rz(0.49219254) q[0];
rz(-2.4824712) q[1];
sx q[1];
rz(-0.4159795) q[1];
sx q[1];
rz(0.98512828) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8644442) q[0];
sx q[0];
rz(-1.5168958) q[0];
sx q[0];
rz(-0.076907579) q[0];
x q[1];
rz(1.6798626) q[2];
sx q[2];
rz(-1.1972053) q[2];
sx q[2];
rz(-1.9337758) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0790365) q[1];
sx q[1];
rz(-0.46508671) q[1];
sx q[1];
rz(-2.5885838) q[1];
rz(-pi) q[2];
rz(1.5268232) q[3];
sx q[3];
rz(-0.73426651) q[3];
sx q[3];
rz(1.6787613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.27348406) q[2];
sx q[2];
rz(-2.2874139) q[2];
sx q[2];
rz(-1.2218529) q[2];
rz(-2.6814804) q[3];
sx q[3];
rz(-1.2686138) q[3];
sx q[3];
rz(2.4251895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7465794) q[0];
sx q[0];
rz(-1.6759251) q[0];
sx q[0];
rz(1.2431086) q[0];
rz(1.5383491) q[1];
sx q[1];
rz(-1.0335045) q[1];
sx q[1];
rz(-0.43362591) q[1];
rz(1.4547841) q[2];
sx q[2];
rz(-0.83509904) q[2];
sx q[2];
rz(-0.72088045) q[2];
rz(-3.1258413) q[3];
sx q[3];
rz(-1.2917662) q[3];
sx q[3];
rz(0.36242604) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
