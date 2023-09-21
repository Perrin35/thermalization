OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.38874415) q[0];
sx q[0];
rz(3.677877) q[0];
sx q[0];
rz(10.372547) q[0];
rz(1.8127958) q[1];
sx q[1];
rz(-1.2674018) q[1];
sx q[1];
rz(-1.0277494) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42067895) q[0];
sx q[0];
rz(-1.5936216) q[0];
sx q[0];
rz(2.7657397) q[0];
x q[1];
rz(-2.920354) q[2];
sx q[2];
rz(-1.0916296) q[2];
sx q[2];
rz(2.0146807) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6537522) q[1];
sx q[1];
rz(-2.2334705) q[1];
sx q[1];
rz(-2.0648048) q[1];
rz(-pi) q[2];
rz(0.58524744) q[3];
sx q[3];
rz(-0.71551502) q[3];
sx q[3];
rz(-1.9714718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1296922) q[2];
sx q[2];
rz(-1.4346069) q[2];
sx q[2];
rz(2.091308) q[2];
rz(2.0283279) q[3];
sx q[3];
rz(-2.2498825) q[3];
sx q[3];
rz(0.068107001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.072409078) q[0];
sx q[0];
rz(-1.2453112) q[0];
sx q[0];
rz(-2.8438399) q[0];
rz(0.61966664) q[1];
sx q[1];
rz(-2.1344118) q[1];
sx q[1];
rz(2.0334977) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8544234) q[0];
sx q[0];
rz(-0.65432917) q[0];
sx q[0];
rz(-1.2229344) q[0];
x q[1];
rz(-0.96252243) q[2];
sx q[2];
rz(-1.2108742) q[2];
sx q[2];
rz(-2.4134709) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.88156466) q[1];
sx q[1];
rz(-1.9238872) q[1];
sx q[1];
rz(-0.41207037) q[1];
x q[2];
rz(-2.8262029) q[3];
sx q[3];
rz(-2.4139801) q[3];
sx q[3];
rz(1.7244347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6790598) q[2];
sx q[2];
rz(-2.1647537) q[2];
sx q[2];
rz(-2.1662946) q[2];
rz(0.9179999) q[3];
sx q[3];
rz(-1.2851597) q[3];
sx q[3];
rz(0.29618922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.179203) q[0];
sx q[0];
rz(-0.85150349) q[0];
sx q[0];
rz(-2.5986824) q[0];
rz(-2.2593598) q[1];
sx q[1];
rz(-2.0062607) q[1];
sx q[1];
rz(-0.96484819) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1436413) q[0];
sx q[0];
rz(-0.43295857) q[0];
sx q[0];
rz(-0.71822449) q[0];
rz(-pi) q[1];
rz(0.29693551) q[2];
sx q[2];
rz(-1.529971) q[2];
sx q[2];
rz(-2.9253935) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.45005349) q[1];
sx q[1];
rz(-1.69194) q[1];
sx q[1];
rz(0.064784592) q[1];
rz(-pi) q[2];
rz(-2.9433555) q[3];
sx q[3];
rz(-1.9803515) q[3];
sx q[3];
rz(0.73892456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0470011) q[2];
sx q[2];
rz(-2.5307405) q[2];
sx q[2];
rz(1.1331406) q[2];
rz(0.23162332) q[3];
sx q[3];
rz(-1.2730205) q[3];
sx q[3];
rz(2.384322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.27947458) q[0];
sx q[0];
rz(-0.010443895) q[0];
sx q[0];
rz(1.7650771) q[0];
rz(2.6230985) q[1];
sx q[1];
rz(-1.8771749) q[1];
sx q[1];
rz(2.8994697) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74894416) q[0];
sx q[0];
rz(-0.90396515) q[0];
sx q[0];
rz(1.3632266) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5838727) q[2];
sx q[2];
rz(-1.1411975) q[2];
sx q[2];
rz(2.122288) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.33430797) q[1];
sx q[1];
rz(-1.1242928) q[1];
sx q[1];
rz(1.7590894) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4297156) q[3];
sx q[3];
rz(-1.1564848) q[3];
sx q[3];
rz(0.29841081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.018192856) q[2];
sx q[2];
rz(-0.96863666) q[2];
sx q[2];
rz(3.0467395) q[2];
rz(1.799396) q[3];
sx q[3];
rz(-1.7443402) q[3];
sx q[3];
rz(-0.43911394) q[3];
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
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78385329) q[0];
sx q[0];
rz(-1.8308324) q[0];
sx q[0];
rz(2.7807996) q[0];
rz(1.3882673) q[1];
sx q[1];
rz(-1.3307945) q[1];
sx q[1];
rz(-1.1345908) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15123385) q[0];
sx q[0];
rz(-0.84999527) q[0];
sx q[0];
rz(-0.6806586) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1410494) q[2];
sx q[2];
rz(-1.2417214) q[2];
sx q[2];
rz(-0.82440257) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1281631) q[1];
sx q[1];
rz(-1.1268106) q[1];
sx q[1];
rz(1.7432937) q[1];
rz(-0.43679045) q[3];
sx q[3];
rz(-1.4947287) q[3];
sx q[3];
rz(1.2341183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1363042) q[2];
sx q[2];
rz(-1.4207999) q[2];
sx q[2];
rz(-2.5689382) q[2];
rz(-0.92875656) q[3];
sx q[3];
rz(-0.52162617) q[3];
sx q[3];
rz(1.9740392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3271493) q[0];
sx q[0];
rz(-1.0961908) q[0];
sx q[0];
rz(-1.8922528) q[0];
rz(-1.918474) q[1];
sx q[1];
rz(-1.5253116) q[1];
sx q[1];
rz(1.1522326) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0056861176) q[0];
sx q[0];
rz(-1.5581157) q[0];
sx q[0];
rz(-1.6760875) q[0];
rz(-1.9828898) q[2];
sx q[2];
rz(-2.0372143) q[2];
sx q[2];
rz(2.7127624) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9705379) q[1];
sx q[1];
rz(-1.4661403) q[1];
sx q[1];
rz(-2.9403951) q[1];
rz(-0.67752083) q[3];
sx q[3];
rz(-2.065425) q[3];
sx q[3];
rz(2.236931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46889177) q[2];
sx q[2];
rz(-1.2142618) q[2];
sx q[2];
rz(0.99096283) q[2];
rz(0.64783603) q[3];
sx q[3];
rz(-0.9655374) q[3];
sx q[3];
rz(-2.794054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8836442) q[0];
sx q[0];
rz(-2.914496) q[0];
sx q[0];
rz(0.062285475) q[0];
rz(2.9557872) q[1];
sx q[1];
rz(-1.6848911) q[1];
sx q[1];
rz(-2.7468162) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61790723) q[0];
sx q[0];
rz(-0.46645188) q[0];
sx q[0];
rz(2.594069) q[0];
rz(-pi) q[1];
rz(-1.8918745) q[2];
sx q[2];
rz(-1.9294538) q[2];
sx q[2];
rz(-1.1527588) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3697606) q[1];
sx q[1];
rz(-2.6409915) q[1];
sx q[1];
rz(-2.6066149) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9881334) q[3];
sx q[3];
rz(-1.6073265) q[3];
sx q[3];
rz(-0.032144459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.64849598) q[2];
sx q[2];
rz(-2.811537) q[2];
sx q[2];
rz(-0.27302343) q[2];
rz(1.8388883) q[3];
sx q[3];
rz(-1.3132934) q[3];
sx q[3];
rz(-0.31204143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39772314) q[0];
sx q[0];
rz(-2.0070772) q[0];
sx q[0];
rz(1.2851108) q[0];
rz(1.5015191) q[1];
sx q[1];
rz(-1.7506426) q[1];
sx q[1];
rz(1.8008908) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9467981) q[0];
sx q[0];
rz(-1.9614944) q[0];
sx q[0];
rz(-2.7892053) q[0];
rz(3.0807207) q[2];
sx q[2];
rz(-0.69673046) q[2];
sx q[2];
rz(2.3827162) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.33412877) q[1];
sx q[1];
rz(-2.6263116) q[1];
sx q[1];
rz(0.79828243) q[1];
rz(-pi) q[2];
rz(0.96745456) q[3];
sx q[3];
rz(-0.52166044) q[3];
sx q[3];
rz(1.1286917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1061219) q[2];
sx q[2];
rz(-0.80646986) q[2];
sx q[2];
rz(-2.1179874) q[2];
rz(-2.9566531) q[3];
sx q[3];
rz(-2.7513294) q[3];
sx q[3];
rz(0.24188724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0446562) q[0];
sx q[0];
rz(-0.99656314) q[0];
sx q[0];
rz(1.6212844) q[0];
rz(-2.8114491) q[1];
sx q[1];
rz(-1.9338927) q[1];
sx q[1];
rz(-2.3044589) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9159106) q[0];
sx q[0];
rz(-1.3689694) q[0];
sx q[0];
rz(0.68737824) q[0];
rz(-pi) q[1];
rz(-0.57640055) q[2];
sx q[2];
rz(-1.8216368) q[2];
sx q[2];
rz(-2.5121411) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3763807) q[1];
sx q[1];
rz(-1.7227255) q[1];
sx q[1];
rz(3.0357749) q[1];
rz(-0.29104851) q[3];
sx q[3];
rz(-1.0927199) q[3];
sx q[3];
rz(0.92791286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.60124406) q[2];
sx q[2];
rz(-1.0849755) q[2];
sx q[2];
rz(-2.9619651) q[2];
rz(-2.1458697) q[3];
sx q[3];
rz(-1.2487753) q[3];
sx q[3];
rz(-1.3109591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3778465) q[0];
sx q[0];
rz(-2.7959931) q[0];
sx q[0];
rz(1.0572222) q[0];
rz(3.0341042) q[1];
sx q[1];
rz(-1.8881533) q[1];
sx q[1];
rz(0.9799788) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5883334) q[0];
sx q[0];
rz(-0.13561121) q[0];
sx q[0];
rz(1.1853663) q[0];
x q[1];
rz(-0.76358958) q[2];
sx q[2];
rz(-0.66765235) q[2];
sx q[2];
rz(-1.9612567) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.53303888) q[1];
sx q[1];
rz(-1.8095008) q[1];
sx q[1];
rz(2.0134316) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8944593) q[3];
sx q[3];
rz(-1.4753046) q[3];
sx q[3];
rz(3.132706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3048627) q[2];
sx q[2];
rz(-1.9367846) q[2];
sx q[2];
rz(-1.0277964) q[2];
rz(1.7547539) q[3];
sx q[3];
rz(-1.3091062) q[3];
sx q[3];
rz(0.28361472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86823157) q[0];
sx q[0];
rz(-1.0366476) q[0];
sx q[0];
rz(-1.114053) q[0];
rz(2.3095619) q[1];
sx q[1];
rz(-2.6770626) q[1];
sx q[1];
rz(-2.4774036) q[1];
rz(-1.8503415) q[2];
sx q[2];
rz(-2.5646979) q[2];
sx q[2];
rz(-0.23451351) q[2];
rz(-2.5084393) q[3];
sx q[3];
rz(-0.59637759) q[3];
sx q[3];
rz(-0.37806088) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];