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
rz(1.8019851) q[0];
sx q[0];
rz(-2.7924502) q[0];
sx q[0];
rz(-2.1139297) q[0];
rz(-0.034962058) q[1];
sx q[1];
rz(5.265994) q[1];
sx q[1];
rz(7.979402) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48533121) q[0];
sx q[0];
rz(-2.4527271) q[0];
sx q[0];
rz(2.5883677) q[0];
rz(-0.18741636) q[2];
sx q[2];
rz(-0.18067154) q[2];
sx q[2];
rz(-2.540321) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6830292) q[1];
sx q[1];
rz(-1.7818319) q[1];
sx q[1];
rz(-0.33331916) q[1];
x q[2];
rz(-1.2054005) q[3];
sx q[3];
rz(-0.8626079) q[3];
sx q[3];
rz(-2.303249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.33862904) q[2];
sx q[2];
rz(-0.45946071) q[2];
sx q[2];
rz(0.53172338) q[2];
rz(-1.7470597) q[3];
sx q[3];
rz(-1.8683542) q[3];
sx q[3];
rz(1.9664221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26674536) q[0];
sx q[0];
rz(-1.5044455) q[0];
sx q[0];
rz(-2.3496085) q[0];
rz(0.92164552) q[1];
sx q[1];
rz(-1.6519203) q[1];
sx q[1];
rz(-2.0335061) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9757248) q[0];
sx q[0];
rz(-2.168403) q[0];
sx q[0];
rz(0.13522603) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5032866) q[2];
sx q[2];
rz(-0.6246399) q[2];
sx q[2];
rz(-0.72062525) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0217359) q[1];
sx q[1];
rz(-0.83470063) q[1];
sx q[1];
rz(2.0935007) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4008303) q[3];
sx q[3];
rz(-2.8208591) q[3];
sx q[3];
rz(0.50839822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.38829982) q[2];
sx q[2];
rz(-2.1123977) q[2];
sx q[2];
rz(-1.252906) q[2];
rz(-2.5168354) q[3];
sx q[3];
rz(-1.2478991) q[3];
sx q[3];
rz(-2.6784082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4598292) q[0];
sx q[0];
rz(-2.9637931) q[0];
sx q[0];
rz(2.8681927) q[0];
rz(-3.0699442) q[1];
sx q[1];
rz(-2.007808) q[1];
sx q[1];
rz(-1.6260446) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0833383) q[0];
sx q[0];
rz(-1.951955) q[0];
sx q[0];
rz(2.5883771) q[0];
rz(2.6726682) q[2];
sx q[2];
rz(-1.8919626) q[2];
sx q[2];
rz(1.4844446) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3699941) q[1];
sx q[1];
rz(-0.99111667) q[1];
sx q[1];
rz(-1.3808151) q[1];
x q[2];
rz(0.97730009) q[3];
sx q[3];
rz(-1.3588828) q[3];
sx q[3];
rz(-0.43928543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.36149725) q[2];
sx q[2];
rz(-2.4269035) q[2];
sx q[2];
rz(-1.7717465) q[2];
rz(2.0731549) q[3];
sx q[3];
rz(-1.3911824) q[3];
sx q[3];
rz(0.95503241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2925828) q[0];
sx q[0];
rz(-2.7136901) q[0];
sx q[0];
rz(0.12246116) q[0];
rz(3.124584) q[1];
sx q[1];
rz(-1.6288792) q[1];
sx q[1];
rz(0.88517991) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95851446) q[0];
sx q[0];
rz(-1.4349271) q[0];
sx q[0];
rz(-2.0625173) q[0];
rz(-pi) q[1];
x q[1];
rz(0.34220747) q[2];
sx q[2];
rz(-2.2762083) q[2];
sx q[2];
rz(-0.17555412) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.95751132) q[1];
sx q[1];
rz(-2.5169417) q[1];
sx q[1];
rz(0.84926012) q[1];
x q[2];
rz(1.3766798) q[3];
sx q[3];
rz(-0.87820617) q[3];
sx q[3];
rz(-0.10759456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7077606) q[2];
sx q[2];
rz(-1.8261352) q[2];
sx q[2];
rz(-0.07746499) q[2];
rz(2.7205983) q[3];
sx q[3];
rz(-1.1869895) q[3];
sx q[3];
rz(-0.25119701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0780599) q[0];
sx q[0];
rz(-3.1227626) q[0];
sx q[0];
rz(-0.68191648) q[0];
rz(-2.6497427) q[1];
sx q[1];
rz(-1.0268772) q[1];
sx q[1];
rz(1.9815725) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9043961) q[0];
sx q[0];
rz(-0.86441308) q[0];
sx q[0];
rz(-0.24817962) q[0];
rz(1.5142646) q[2];
sx q[2];
rz(-1.0011359) q[2];
sx q[2];
rz(0.73318549) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8012538) q[1];
sx q[1];
rz(-1.7427708) q[1];
sx q[1];
rz(0.94229631) q[1];
rz(-pi) q[2];
rz(-2.8468708) q[3];
sx q[3];
rz(-1.975276) q[3];
sx q[3];
rz(1.3724788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.457095) q[2];
sx q[2];
rz(-0.86898494) q[2];
sx q[2];
rz(1.0082461) q[2];
rz(-1.6393939) q[3];
sx q[3];
rz(-2.0366171) q[3];
sx q[3];
rz(-2.8777299) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.098792583) q[0];
sx q[0];
rz(-1.5941987) q[0];
sx q[0];
rz(-2.7368326) q[0];
rz(2.3233991) q[1];
sx q[1];
rz(-1.300756) q[1];
sx q[1];
rz(2.4249446) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7966864) q[0];
sx q[0];
rz(-0.27936882) q[0];
sx q[0];
rz(-0.035682364) q[0];
rz(-2.0990853) q[2];
sx q[2];
rz(-2.2620438) q[2];
sx q[2];
rz(1.7886358) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4617052) q[1];
sx q[1];
rz(-2.6033834) q[1];
sx q[1];
rz(2.6160282) q[1];
rz(-pi) q[2];
rz(-0.53933177) q[3];
sx q[3];
rz(-1.2507273) q[3];
sx q[3];
rz(2.0254997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3396259) q[2];
sx q[2];
rz(-2.5012987) q[2];
sx q[2];
rz(-2.4862945) q[2];
rz(-0.35704923) q[3];
sx q[3];
rz(-1.7506295) q[3];
sx q[3];
rz(-2.6220139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94012564) q[0];
sx q[0];
rz(-0.31620142) q[0];
sx q[0];
rz(2.1215718) q[0];
rz(2.5992498) q[1];
sx q[1];
rz(-2.0261363) q[1];
sx q[1];
rz(-1.7592336) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6665628) q[0];
sx q[0];
rz(-2.2336166) q[0];
sx q[0];
rz(0.7866089) q[0];
rz(2.9303161) q[2];
sx q[2];
rz(-1.7045492) q[2];
sx q[2];
rz(1.1115896) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0918563) q[1];
sx q[1];
rz(-2.0309134) q[1];
sx q[1];
rz(3.0541944) q[1];
rz(-0.027023836) q[3];
sx q[3];
rz(-3.0579429) q[3];
sx q[3];
rz(0.25318709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.072784) q[2];
sx q[2];
rz(-1.5280318) q[2];
sx q[2];
rz(2.7911348) q[2];
rz(-2.6751878) q[3];
sx q[3];
rz(-2.2240708) q[3];
sx q[3];
rz(2.1980227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7826409) q[0];
sx q[0];
rz(-2.751001) q[0];
sx q[0];
rz(2.1851831) q[0];
rz(-1.4495173) q[1];
sx q[1];
rz(-2.3608975) q[1];
sx q[1];
rz(0.67109674) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3365375) q[0];
sx q[0];
rz(-1.4292246) q[0];
sx q[0];
rz(0.086351589) q[0];
x q[1];
rz(2.5726357) q[2];
sx q[2];
rz(-1.3929741) q[2];
sx q[2];
rz(-0.59019719) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.073428) q[1];
sx q[1];
rz(-2.7824017) q[1];
sx q[1];
rz(1.3884307) q[1];
rz(-0.70399474) q[3];
sx q[3];
rz(-1.7592832) q[3];
sx q[3];
rz(-0.89787591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.87390071) q[2];
sx q[2];
rz(-2.716422) q[2];
sx q[2];
rz(1.1262061) q[2];
rz(-2.8113484) q[3];
sx q[3];
rz(-1.4010022) q[3];
sx q[3];
rz(-0.11708524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.0983593) q[0];
sx q[0];
rz(-0.57906228) q[0];
sx q[0];
rz(-0.057057127) q[0];
rz(0.13380274) q[1];
sx q[1];
rz(-1.0452784) q[1];
sx q[1];
rz(-0.71279508) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.063449115) q[0];
sx q[0];
rz(-0.3219372) q[0];
sx q[0];
rz(2.9418895) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1634689) q[2];
sx q[2];
rz(-1.9453269) q[2];
sx q[2];
rz(-0.56331149) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7851398) q[1];
sx q[1];
rz(-0.63842652) q[1];
sx q[1];
rz(0.70968763) q[1];
rz(-pi) q[2];
rz(3.0607575) q[3];
sx q[3];
rz(-2.2634215) q[3];
sx q[3];
rz(2.5580542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.628525) q[2];
sx q[2];
rz(-1.8551989) q[2];
sx q[2];
rz(3.0600424) q[2];
rz(-1.2201803) q[3];
sx q[3];
rz(-2.8704075) q[3];
sx q[3];
rz(-1.0903821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9995025) q[0];
sx q[0];
rz(-1.6195848) q[0];
sx q[0];
rz(1.5600486) q[0];
rz(2.0060495) q[1];
sx q[1];
rz(-1.138843) q[1];
sx q[1];
rz(-1.9948237) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97871507) q[0];
sx q[0];
rz(-0.93746569) q[0];
sx q[0];
rz(0.94453728) q[0];
rz(-pi) q[1];
rz(0.79939553) q[2];
sx q[2];
rz(-2.172432) q[2];
sx q[2];
rz(-0.55527273) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9830977) q[1];
sx q[1];
rz(-1.4870411) q[1];
sx q[1];
rz(-0.2106481) q[1];
x q[2];
rz(0.087298079) q[3];
sx q[3];
rz(-1.0074769) q[3];
sx q[3];
rz(-0.47730744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.59209383) q[2];
sx q[2];
rz(-1.8213976) q[2];
sx q[2];
rz(-0.36592323) q[2];
rz(-2.7538815) q[3];
sx q[3];
rz(-1.1185027) q[3];
sx q[3];
rz(1.4661192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53030071) q[0];
sx q[0];
rz(-2.3983751) q[0];
sx q[0];
rz(-0.30839738) q[0];
rz(-2.3920234) q[1];
sx q[1];
rz(-2.0105965) q[1];
sx q[1];
rz(-1.2731193) q[1];
rz(-2.0012435) q[2];
sx q[2];
rz(-0.7067718) q[2];
sx q[2];
rz(1.5511647) q[2];
rz(-0.6324296) q[3];
sx q[3];
rz(-1.2865744) q[3];
sx q[3];
rz(-2.1435973) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
