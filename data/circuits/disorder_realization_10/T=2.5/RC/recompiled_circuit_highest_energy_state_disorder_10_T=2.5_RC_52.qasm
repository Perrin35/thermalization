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
rz(-0.69775692) q[0];
sx q[0];
rz(-1.193576) q[0];
sx q[0];
rz(2.9370263) q[0];
rz(-0.76072955) q[1];
sx q[1];
rz(-1.3890356) q[1];
sx q[1];
rz(1.3995481) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22628173) q[0];
sx q[0];
rz(-1.4797512) q[0];
sx q[0];
rz(1.7928726) q[0];
x q[1];
rz(-0.40488196) q[2];
sx q[2];
rz(-1.153115) q[2];
sx q[2];
rz(1.6987125) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.6217566) q[1];
sx q[1];
rz(-2.7108828) q[1];
sx q[1];
rz(2.9309896) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8329323) q[3];
sx q[3];
rz(-0.7334358) q[3];
sx q[3];
rz(1.6876458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.43510398) q[2];
sx q[2];
rz(-0.032328345) q[2];
sx q[2];
rz(0.48933634) q[2];
rz(1.0232183) q[3];
sx q[3];
rz(-3.1232941) q[3];
sx q[3];
rz(-2.0160915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.045687549) q[0];
sx q[0];
rz(-0.65653312) q[0];
sx q[0];
rz(-0.80279654) q[0];
rz(-0.071391694) q[1];
sx q[1];
rz(-2.8749021) q[1];
sx q[1];
rz(-3.0838222) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36299713) q[0];
sx q[0];
rz(-2.5311806) q[0];
sx q[0];
rz(1.1108062) q[0];
rz(2.8550496) q[2];
sx q[2];
rz(-1.4689969) q[2];
sx q[2];
rz(1.2004146) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3398185) q[1];
sx q[1];
rz(-1.8027824) q[1];
sx q[1];
rz(1.07527) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14795078) q[3];
sx q[3];
rz(-0.86455621) q[3];
sx q[3];
rz(1.0195635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.94566655) q[2];
sx q[2];
rz(-1.095093) q[2];
sx q[2];
rz(-1.2954953) q[2];
rz(0.96674353) q[3];
sx q[3];
rz(-0.77015489) q[3];
sx q[3];
rz(-2.3841592) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029723786) q[0];
sx q[0];
rz(-1.7787378) q[0];
sx q[0];
rz(1.7080074) q[0];
rz(-0.068610527) q[1];
sx q[1];
rz(-1.5741293) q[1];
sx q[1];
rz(0.57919085) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60459701) q[0];
sx q[0];
rz(-1.662435) q[0];
sx q[0];
rz(3.0834404) q[0];
rz(-pi) q[1];
rz(0.0036598031) q[2];
sx q[2];
rz(-0.47572593) q[2];
sx q[2];
rz(-1.9329247) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3855558) q[1];
sx q[1];
rz(-1.7969404) q[1];
sx q[1];
rz(-1.6075587) q[1];
rz(-pi) q[2];
rz(-2.0491568) q[3];
sx q[3];
rz(-1.98171) q[3];
sx q[3];
rz(-2.876407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8706943) q[2];
sx q[2];
rz(-1.8879994) q[2];
sx q[2];
rz(-2.9581621) q[2];
rz(2.3373248) q[3];
sx q[3];
rz(-0.99884123) q[3];
sx q[3];
rz(-2.6075294) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19876984) q[0];
sx q[0];
rz(-3.0267921) q[0];
sx q[0];
rz(-0.59952366) q[0];
rz(2.7291258) q[1];
sx q[1];
rz(-0.020655276) q[1];
sx q[1];
rz(-2.1309158) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0507999) q[0];
sx q[0];
rz(-0.37380344) q[0];
sx q[0];
rz(0.35275491) q[0];
rz(-1.640219) q[2];
sx q[2];
rz(-0.63328082) q[2];
sx q[2];
rz(-0.4236003) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7059203) q[1];
sx q[1];
rz(-1.3002035) q[1];
sx q[1];
rz(2.1822004) q[1];
x q[2];
rz(-1.7092199) q[3];
sx q[3];
rz(-2.3218621) q[3];
sx q[3];
rz(-2.8734796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7708873) q[2];
sx q[2];
rz(-0.34660307) q[2];
sx q[2];
rz(2.8240805) q[2];
rz(-2.5192449) q[3];
sx q[3];
rz(-2.205866) q[3];
sx q[3];
rz(0.58421016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36251003) q[0];
sx q[0];
rz(-0.91479397) q[0];
sx q[0];
rz(1.0269748) q[0];
rz(2.5912071) q[1];
sx q[1];
rz(-3.0774979) q[1];
sx q[1];
rz(1.9245573) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7893697) q[0];
sx q[0];
rz(-0.93607157) q[0];
sx q[0];
rz(2.288649) q[0];
x q[1];
rz(-3.0116357) q[2];
sx q[2];
rz(-0.54207793) q[2];
sx q[2];
rz(0.96058577) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.74351013) q[1];
sx q[1];
rz(-0.71894246) q[1];
sx q[1];
rz(1.8330857) q[1];
rz(-pi) q[2];
rz(0.86037614) q[3];
sx q[3];
rz(-1.2902707) q[3];
sx q[3];
rz(-0.70395148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7050742) q[2];
sx q[2];
rz(-1.3631835) q[2];
sx q[2];
rz(-0.77318937) q[2];
rz(0.1117205) q[3];
sx q[3];
rz(-1.819928) q[3];
sx q[3];
rz(-0.93938655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47950995) q[0];
sx q[0];
rz(-0.22308068) q[0];
sx q[0];
rz(-2.7146085) q[0];
rz(2.2110979) q[1];
sx q[1];
rz(-0.016914802) q[1];
sx q[1];
rz(-2.6771136) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8278566) q[0];
sx q[0];
rz(-1.858874) q[0];
sx q[0];
rz(1.366607) q[0];
rz(-pi) q[1];
rz(-0.97067483) q[2];
sx q[2];
rz(-1.8541186) q[2];
sx q[2];
rz(1.4221869) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0911922) q[1];
sx q[1];
rz(-2.770863) q[1];
sx q[1];
rz(-0.099845366) q[1];
x q[2];
rz(-1.1464981) q[3];
sx q[3];
rz(-1.5992377) q[3];
sx q[3];
rz(-0.44876307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.53133196) q[2];
sx q[2];
rz(-1.8208296) q[2];
sx q[2];
rz(-2.8533234) q[2];
rz(1.0432976) q[3];
sx q[3];
rz(-2.5259924) q[3];
sx q[3];
rz(-2.402795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4814602) q[0];
sx q[0];
rz(-1.2876502) q[0];
sx q[0];
rz(-0.75755358) q[0];
rz(-0.07269147) q[1];
sx q[1];
rz(-0.025765954) q[1];
sx q[1];
rz(0.048197897) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2295727) q[0];
sx q[0];
rz(-0.98668146) q[0];
sx q[0];
rz(0.035365625) q[0];
rz(-pi) q[1];
rz(-1.2343302) q[2];
sx q[2];
rz(-2.5162553) q[2];
sx q[2];
rz(2.4246374) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9513598) q[1];
sx q[1];
rz(-2.1589365) q[1];
sx q[1];
rz(-2.6201671) q[1];
rz(-pi) q[2];
rz(1.1009351) q[3];
sx q[3];
rz(-1.5239232) q[3];
sx q[3];
rz(2.8702503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4797719) q[2];
sx q[2];
rz(-1.6419819) q[2];
sx q[2];
rz(-3.0721967) q[2];
rz(1.5540468) q[3];
sx q[3];
rz(-2.3543251) q[3];
sx q[3];
rz(0.21849304) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3875535) q[0];
sx q[0];
rz(-1.0559005) q[0];
sx q[0];
rz(-1.4132502) q[0];
rz(-2.3130401) q[1];
sx q[1];
rz(-3.0998402) q[1];
sx q[1];
rz(-0.54263306) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6814282) q[0];
sx q[0];
rz(-1.6896473) q[0];
sx q[0];
rz(-2.9959841) q[0];
rz(1.9285525) q[2];
sx q[2];
rz(-1.2160436) q[2];
sx q[2];
rz(2.8194129) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.97510135) q[1];
sx q[1];
rz(-1.5169739) q[1];
sx q[1];
rz(0.11589284) q[1];
rz(1.1922791) q[3];
sx q[3];
rz(-1.9782441) q[3];
sx q[3];
rz(-2.440314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1594306) q[2];
sx q[2];
rz(-2.7320778) q[2];
sx q[2];
rz(2.8816667) q[2];
rz(-2.121117) q[3];
sx q[3];
rz(-0.25769886) q[3];
sx q[3];
rz(2.3475588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6771773) q[0];
sx q[0];
rz(-0.18615119) q[0];
sx q[0];
rz(1.6625241) q[0];
rz(1.6089449) q[1];
sx q[1];
rz(-2.1023991) q[1];
sx q[1];
rz(-2.398568) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0940762) q[0];
sx q[0];
rz(-2.0819944) q[0];
sx q[0];
rz(0.1316977) q[0];
rz(-2.7588927) q[2];
sx q[2];
rz(-2.0819132) q[2];
sx q[2];
rz(1.9267043) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8919395) q[1];
sx q[1];
rz(-1.5158476) q[1];
sx q[1];
rz(-1.5049388) q[1];
x q[2];
rz(2.8554932) q[3];
sx q[3];
rz(-0.10933441) q[3];
sx q[3];
rz(-0.63150333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.6741901) q[2];
sx q[2];
rz(-2.3154066) q[2];
sx q[2];
rz(0.80545938) q[2];
rz(1.6921267) q[3];
sx q[3];
rz(-1.9129246) q[3];
sx q[3];
rz(0.8031556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8886803) q[0];
sx q[0];
rz(-0.54179931) q[0];
sx q[0];
rz(2.3895277) q[0];
rz(-1.9883142) q[1];
sx q[1];
rz(-0.88429943) q[1];
sx q[1];
rz(0.28336743) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7572989) q[0];
sx q[0];
rz(-1.6597676) q[0];
sx q[0];
rz(0.53953895) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7287909) q[2];
sx q[2];
rz(-2.3834977) q[2];
sx q[2];
rz(1.5811063) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1438469) q[1];
sx q[1];
rz(-2.4334868) q[1];
sx q[1];
rz(-2.6216402) q[1];
rz(-pi) q[2];
x q[2];
rz(2.649077) q[3];
sx q[3];
rz(-2.4253143) q[3];
sx q[3];
rz(-0.034958358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.89455426) q[2];
sx q[2];
rz(-3.0587695) q[2];
sx q[2];
rz(1.7029597) q[2];
rz(0.29397193) q[3];
sx q[3];
rz(-0.014466244) q[3];
sx q[3];
rz(-1.0283874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1080078) q[0];
sx q[0];
rz(-1.7243732) q[0];
sx q[0];
rz(1.617817) q[0];
rz(2.602018) q[1];
sx q[1];
rz(-2.3537666) q[1];
sx q[1];
rz(-2.981577) q[1];
rz(2.9931184) q[2];
sx q[2];
rz(-1.21429) q[2];
sx q[2];
rz(0.16333632) q[2];
rz(-2.5359375) q[3];
sx q[3];
rz(-0.45481053) q[3];
sx q[3];
rz(-0.74301471) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
