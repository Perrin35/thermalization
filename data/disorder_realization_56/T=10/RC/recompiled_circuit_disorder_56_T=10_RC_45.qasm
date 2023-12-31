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
rz(2.1138432) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2078903) q[0];
sx q[0];
rz(-2.7650802) q[0];
sx q[0];
rz(3.0794789) q[0];
rz(1.0814632) q[2];
sx q[2];
rz(-1.3748193) q[2];
sx q[2];
rz(-2.8010362) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9103968) q[1];
sx q[1];
rz(-0.80363552) q[1];
sx q[1];
rz(0.5459783) q[1];
x q[2];
rz(2.5563452) q[3];
sx q[3];
rz(-2.4260776) q[3];
sx q[3];
rz(1.1701208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1296922) q[2];
sx q[2];
rz(-1.4346069) q[2];
sx q[2];
rz(1.0502846) q[2];
rz(2.0283279) q[3];
sx q[3];
rz(-0.89171019) q[3];
sx q[3];
rz(-0.068107001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.072409078) q[0];
sx q[0];
rz(-1.8962815) q[0];
sx q[0];
rz(-0.29775277) q[0];
rz(-2.521926) q[1];
sx q[1];
rz(-2.1344118) q[1];
sx q[1];
rz(2.0334977) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8544234) q[0];
sx q[0];
rz(-0.65432917) q[0];
sx q[0];
rz(-1.9186583) q[0];
rz(-pi) q[1];
rz(-0.98845311) q[2];
sx q[2];
rz(-0.69500143) q[2];
sx q[2];
rz(-0.37441355) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.83924343) q[1];
sx q[1];
rz(-1.1855372) q[1];
sx q[1];
rz(-1.1883931) q[1];
rz(-0.31538972) q[3];
sx q[3];
rz(-0.72761256) q[3];
sx q[3];
rz(1.7244347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6790598) q[2];
sx q[2];
rz(-0.97683895) q[2];
sx q[2];
rz(0.97529808) q[2];
rz(-2.2235928) q[3];
sx q[3];
rz(-1.2851597) q[3];
sx q[3];
rz(-2.8454034) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.179203) q[0];
sx q[0];
rz(-0.85150349) q[0];
sx q[0];
rz(-0.54291022) q[0];
rz(0.88223282) q[1];
sx q[1];
rz(-2.0062607) q[1];
sx q[1];
rz(2.1767445) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99795139) q[0];
sx q[0];
rz(-0.43295857) q[0];
sx q[0];
rz(0.71822449) q[0];
x q[1];
rz(2.8446571) q[2];
sx q[2];
rz(-1.529971) q[2];
sx q[2];
rz(-0.21619913) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1988586) q[1];
sx q[1];
rz(-3.004289) q[1];
sx q[1];
rz(1.0820461) q[1];
x q[2];
rz(1.996702) q[3];
sx q[3];
rz(-0.45255462) q[3];
sx q[3];
rz(1.2061314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0470011) q[2];
sx q[2];
rz(-0.61085218) q[2];
sx q[2];
rz(2.0084521) q[2];
rz(-2.9099693) q[3];
sx q[3];
rz(-1.2730205) q[3];
sx q[3];
rz(-0.75727063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.8621181) q[0];
sx q[0];
rz(-3.1311488) q[0];
sx q[0];
rz(-1.7650771) q[0];
rz(-2.6230985) q[1];
sx q[1];
rz(-1.2644178) q[1];
sx q[1];
rz(-0.24212295) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3926485) q[0];
sx q[0];
rz(-0.90396515) q[0];
sx q[0];
rz(-1.3632266) q[0];
rz(-0.02853407) q[2];
sx q[2];
rz(-2.7118073) q[2];
sx q[2];
rz(2.1536749) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.33430797) q[1];
sx q[1];
rz(-1.1242928) q[1];
sx q[1];
rz(1.3825033) q[1];
rz(-pi) q[2];
rz(-1.4297156) q[3];
sx q[3];
rz(-1.9851079) q[3];
sx q[3];
rz(0.29841081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1233998) q[2];
sx q[2];
rz(-2.172956) q[2];
sx q[2];
rz(3.0467395) q[2];
rz(1.799396) q[3];
sx q[3];
rz(-1.7443402) q[3];
sx q[3];
rz(2.7024787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.3577394) q[0];
sx q[0];
rz(-1.8308324) q[0];
sx q[0];
rz(0.36079303) q[0];
rz(-1.7533253) q[1];
sx q[1];
rz(-1.8107982) q[1];
sx q[1];
rz(-2.0070019) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2312233) q[0];
sx q[0];
rz(-2.0632422) q[0];
sx q[0];
rz(-2.4173196) q[0];
rz(-pi) q[1];
rz(0.88419948) q[2];
sx q[2];
rz(-2.6066385) q[2];
sx q[2];
rz(1.3605489) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.36793383) q[1];
sx q[1];
rz(-1.415167) q[1];
sx q[1];
rz(-0.44981287) q[1];
x q[2];
rz(-0.17825019) q[3];
sx q[3];
rz(-0.44294391) q[3];
sx q[3];
rz(2.6435542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1363042) q[2];
sx q[2];
rz(-1.7207928) q[2];
sx q[2];
rz(-2.5689382) q[2];
rz(2.2128361) q[3];
sx q[3];
rz(-0.52162617) q[3];
sx q[3];
rz(-1.1675534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3271493) q[0];
sx q[0];
rz(-1.0961908) q[0];
sx q[0];
rz(1.8922528) q[0];
rz(-1.918474) q[1];
sx q[1];
rz(-1.5253116) q[1];
sx q[1];
rz(1.1522326) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5778225) q[0];
sx q[0];
rz(-1.4655136) q[0];
sx q[0];
rz(3.1288414) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4695775) q[2];
sx q[2];
rz(-2.5294371) q[2];
sx q[2];
rz(-1.200058) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9705379) q[1];
sx q[1];
rz(-1.6754524) q[1];
sx q[1];
rz(-2.9403951) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96529393) q[3];
sx q[3];
rz(-2.1552342) q[3];
sx q[3];
rz(0.30129978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6727009) q[2];
sx q[2];
rz(-1.9273309) q[2];
sx q[2];
rz(0.99096283) q[2];
rz(0.64783603) q[3];
sx q[3];
rz(-2.1760553) q[3];
sx q[3];
rz(2.794054) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8836442) q[0];
sx q[0];
rz(-2.914496) q[0];
sx q[0];
rz(-3.0793072) q[0];
rz(-2.9557872) q[1];
sx q[1];
rz(-1.6848911) q[1];
sx q[1];
rz(-0.3947765) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5236854) q[0];
sx q[0];
rz(-2.6751408) q[0];
sx q[0];
rz(-0.54752366) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4418418) q[2];
sx q[2];
rz(-0.47669461) q[2];
sx q[2];
rz(1.9112019) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.77183206) q[1];
sx q[1];
rz(-0.50060111) q[1];
sx q[1];
rz(-2.6066149) q[1];
rz(1.6077605) q[3];
sx q[3];
rz(-1.4174403) q[3];
sx q[3];
rz(1.6085898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.64849598) q[2];
sx q[2];
rz(-2.811537) q[2];
sx q[2];
rz(-2.8685692) q[2];
rz(1.3027044) q[3];
sx q[3];
rz(-1.8282993) q[3];
sx q[3];
rz(2.8295512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7438695) q[0];
sx q[0];
rz(-1.1345154) q[0];
sx q[0];
rz(1.2851108) q[0];
rz(1.6400736) q[1];
sx q[1];
rz(-1.39095) q[1];
sx q[1];
rz(-1.3407019) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1788951) q[0];
sx q[0];
rz(-2.6216051) q[0];
sx q[0];
rz(0.87332256) q[0];
rz(-pi) q[1];
rz(-1.6216535) q[2];
sx q[2];
rz(-2.2659781) q[2];
sx q[2];
rz(-0.67957544) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9400085) q[1];
sx q[1];
rz(-1.9218947) q[1];
sx q[1];
rz(1.1854978) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0128485) q[3];
sx q[3];
rz(-1.2841409) q[3];
sx q[3];
rz(-3.0451881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0354707) q[2];
sx q[2];
rz(-2.3351228) q[2];
sx q[2];
rz(-2.1179874) q[2];
rz(-2.9566531) q[3];
sx q[3];
rz(-0.39026323) q[3];
sx q[3];
rz(2.8997054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0446562) q[0];
sx q[0];
rz(-0.99656314) q[0];
sx q[0];
rz(-1.6212844) q[0];
rz(-2.8114491) q[1];
sx q[1];
rz(-1.2076999) q[1];
sx q[1];
rz(2.3044589) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9595813) q[0];
sx q[0];
rz(-0.89996979) q[0];
sx q[0];
rz(1.3120033) q[0];
x q[1];
rz(1.8673973) q[2];
sx q[2];
rz(-2.1269848) q[2];
sx q[2];
rz(-1.1013168) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.76444641) q[1];
sx q[1];
rz(-0.18491491) q[1];
sx q[1];
rz(0.96692337) q[1];
x q[2];
rz(1.0749531) q[3];
sx q[3];
rz(-1.3131724) q[3];
sx q[3];
rz(0.50592929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.60124406) q[2];
sx q[2];
rz(-2.0566172) q[2];
sx q[2];
rz(2.9619651) q[2];
rz(-0.99572292) q[3];
sx q[3];
rz(-1.8928173) q[3];
sx q[3];
rz(-1.3109591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3778465) q[0];
sx q[0];
rz(-0.34559956) q[0];
sx q[0];
rz(-1.0572222) q[0];
rz(-0.10748848) q[1];
sx q[1];
rz(-1.8881533) q[1];
sx q[1];
rz(-2.1616139) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9419132) q[0];
sx q[0];
rz(-1.6964039) q[0];
sx q[0];
rz(0.051254) q[0];
rz(-pi) q[1];
rz(0.51771848) q[2];
sx q[2];
rz(-2.0132408) q[2];
sx q[2];
rz(-2.8874318) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.53303888) q[1];
sx q[1];
rz(-1.3320919) q[1];
sx q[1];
rz(-1.1281611) q[1];
rz(-pi) q[2];
rz(-1.8633217) q[3];
sx q[3];
rz(-0.33698002) q[3];
sx q[3];
rz(1.302759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3048627) q[2];
sx q[2];
rz(-1.9367846) q[2];
sx q[2];
rz(2.1137962) q[2];
rz(-1.7547539) q[3];
sx q[3];
rz(-1.3091062) q[3];
sx q[3];
rz(2.8579779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2733611) q[0];
sx q[0];
rz(-2.1049451) q[0];
sx q[0];
rz(2.0275397) q[0];
rz(0.83203075) q[1];
sx q[1];
rz(-0.46453005) q[1];
sx q[1];
rz(0.66418905) q[1];
rz(1.0118532) q[2];
sx q[2];
rz(-1.7218628) q[2];
sx q[2];
rz(-2.0414258) q[2];
rz(2.6408623) q[3];
sx q[3];
rz(-1.9095608) q[3];
sx q[3];
rz(0.6469971) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
