OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.9956545) q[0];
sx q[0];
rz(-0.50322682) q[0];
sx q[0];
rz(-0.72416645) q[0];
rz(-2.5016298) q[1];
sx q[1];
rz(-2.6115186) q[1];
sx q[1];
rz(0.78483265) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9871702) q[0];
sx q[0];
rz(-1.3599456) q[0];
sx q[0];
rz(0.2984557) q[0];
rz(-2.7230524) q[2];
sx q[2];
rz(-1.6633908) q[2];
sx q[2];
rz(1.6469524) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5074871) q[1];
sx q[1];
rz(-1.3480942) q[1];
sx q[1];
rz(2.1102064) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4943487) q[3];
sx q[3];
rz(-0.41562286) q[3];
sx q[3];
rz(1.7474183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.589754) q[2];
sx q[2];
rz(-1.4171615) q[2];
sx q[2];
rz(0.067967728) q[2];
rz(3.0170278) q[3];
sx q[3];
rz(-0.3228651) q[3];
sx q[3];
rz(-1.386806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2215866) q[0];
sx q[0];
rz(-3.0060372) q[0];
sx q[0];
rz(-2.8979229) q[0];
rz(-0.63175732) q[1];
sx q[1];
rz(-1.7383722) q[1];
sx q[1];
rz(-1.3557281) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5403554) q[0];
sx q[0];
rz(-1.8260801) q[0];
sx q[0];
rz(-3.113494) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83265702) q[2];
sx q[2];
rz(-0.50561935) q[2];
sx q[2];
rz(-2.6308699) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.020097453) q[1];
sx q[1];
rz(-1.1250245) q[1];
sx q[1];
rz(-0.14264588) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59216604) q[3];
sx q[3];
rz(-1.6093996) q[3];
sx q[3];
rz(1.6174699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0791066) q[2];
sx q[2];
rz(-0.9920384) q[2];
sx q[2];
rz(0.24965723) q[2];
rz(2.6349973) q[3];
sx q[3];
rz(-1.6258312) q[3];
sx q[3];
rz(2.8095968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(0.24519414) q[0];
sx q[0];
rz(-1.9165374) q[0];
sx q[0];
rz(0.89843345) q[0];
rz(1.8067182) q[1];
sx q[1];
rz(-1.9060262) q[1];
sx q[1];
rz(1.867884) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1184517) q[0];
sx q[0];
rz(-1.3190735) q[0];
sx q[0];
rz(-1.203042) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31295915) q[2];
sx q[2];
rz(-0.3096146) q[2];
sx q[2];
rz(1.6329873) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.18759218) q[1];
sx q[1];
rz(-0.91916579) q[1];
sx q[1];
rz(2.0152394) q[1];
rz(-pi) q[2];
rz(-2.418163) q[3];
sx q[3];
rz(-1.1838786) q[3];
sx q[3];
rz(-0.5667516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0597824) q[2];
sx q[2];
rz(-0.60478294) q[2];
sx q[2];
rz(0.95345062) q[2];
rz(-0.034514286) q[3];
sx q[3];
rz(-0.78648609) q[3];
sx q[3];
rz(2.9147193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2801441) q[0];
sx q[0];
rz(-2.9252958) q[0];
sx q[0];
rz(0.24818534) q[0];
rz(1.0379627) q[1];
sx q[1];
rz(-1.1231517) q[1];
sx q[1];
rz(-3.0674556) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5384597) q[0];
sx q[0];
rz(-2.5939301) q[0];
sx q[0];
rz(-2.0752226) q[0];
rz(-1.9937236) q[2];
sx q[2];
rz(-0.72407702) q[2];
sx q[2];
rz(-2.1507182) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3809507) q[1];
sx q[1];
rz(-1.6886854) q[1];
sx q[1];
rz(1.2537434) q[1];
x q[2];
rz(1.3965963) q[3];
sx q[3];
rz(-1.3341691) q[3];
sx q[3];
rz(3.0803806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2512102) q[2];
sx q[2];
rz(-0.40428287) q[2];
sx q[2];
rz(3.1029491) q[2];
rz(-0.97366992) q[3];
sx q[3];
rz(-0.49574167) q[3];
sx q[3];
rz(2.8715449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-2.6073109) q[0];
sx q[0];
rz(-1.5357635) q[0];
sx q[0];
rz(-1.3624396) q[0];
rz(0.81659395) q[1];
sx q[1];
rz(-1.2885619) q[1];
sx q[1];
rz(1.1626676) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8487932) q[0];
sx q[0];
rz(-1.6822364) q[0];
sx q[0];
rz(0.029990002) q[0];
x q[1];
rz(0.16935279) q[2];
sx q[2];
rz(-1.8893818) q[2];
sx q[2];
rz(-0.40118518) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3711277) q[1];
sx q[1];
rz(-2.7264997) q[1];
sx q[1];
rz(1.0789372) q[1];
rz(-pi) q[2];
rz(0.43443067) q[3];
sx q[3];
rz(-1.0940922) q[3];
sx q[3];
rz(1.3694976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5148619) q[2];
sx q[2];
rz(-2.0613487) q[2];
sx q[2];
rz(0.14222063) q[2];
rz(2.2375315) q[3];
sx q[3];
rz(-1.8217434) q[3];
sx q[3];
rz(2.952125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(3.0267923) q[0];
sx q[0];
rz(-2.8650706) q[0];
sx q[0];
rz(1.4676771) q[0];
rz(2.5698075) q[1];
sx q[1];
rz(-0.3586868) q[1];
sx q[1];
rz(-0.30803672) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7361569) q[0];
sx q[0];
rz(-1.6653367) q[0];
sx q[0];
rz(0.12673881) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.96254827) q[2];
sx q[2];
rz(-2.0934009) q[2];
sx q[2];
rz(-2.7163497) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.088965485) q[1];
sx q[1];
rz(-1.7672156) q[1];
sx q[1];
rz(-0.88167015) q[1];
rz(2.7932348) q[3];
sx q[3];
rz(-1.4847241) q[3];
sx q[3];
rz(1.100988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.77928153) q[2];
sx q[2];
rz(-1.4004536) q[2];
sx q[2];
rz(1.9936838) q[2];
rz(0.71427304) q[3];
sx q[3];
rz(-2.2439984) q[3];
sx q[3];
rz(2.4961297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4797453) q[0];
sx q[0];
rz(-0.84091887) q[0];
sx q[0];
rz(0.1299783) q[0];
rz(0.030844363) q[1];
sx q[1];
rz(-1.2896616) q[1];
sx q[1];
rz(-0.67108363) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86258793) q[0];
sx q[0];
rz(-0.064084856) q[0];
sx q[0];
rz(-0.57871731) q[0];
x q[1];
rz(-1.3865878) q[2];
sx q[2];
rz(-1.5407908) q[2];
sx q[2];
rz(-1.9629994) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.56747251) q[1];
sx q[1];
rz(-1.6735055) q[1];
sx q[1];
rz(1.9436388) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.89550771) q[3];
sx q[3];
rz(-1.607778) q[3];
sx q[3];
rz(-1.0469588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0188296) q[2];
sx q[2];
rz(-1.5315703) q[2];
sx q[2];
rz(2.5637131) q[2];
rz(0.028586483) q[3];
sx q[3];
rz(-1.8609906) q[3];
sx q[3];
rz(-1.2602497) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1639444) q[0];
sx q[0];
rz(-0.90286911) q[0];
sx q[0];
rz(-0.41241616) q[0];
rz(-1.4498129) q[1];
sx q[1];
rz(-1.7990566) q[1];
sx q[1];
rz(1.9746045) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3229423) q[0];
sx q[0];
rz(-2.4837821) q[0];
sx q[0];
rz(-1.2176745) q[0];
rz(-0.86161676) q[2];
sx q[2];
rz(-1.3590727) q[2];
sx q[2];
rz(-2.784563) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.941277) q[1];
sx q[1];
rz(-1.759316) q[1];
sx q[1];
rz(1.539906) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1498296) q[3];
sx q[3];
rz(-2.095788) q[3];
sx q[3];
rz(0.52469745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.940544) q[2];
sx q[2];
rz(-2.2622435) q[2];
sx q[2];
rz(2.9525625) q[2];
rz(-0.14686251) q[3];
sx q[3];
rz(-0.18460128) q[3];
sx q[3];
rz(1.7485025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.015633164) q[0];
sx q[0];
rz(-1.3110315) q[0];
sx q[0];
rz(-0.92700672) q[0];
rz(-1.3828297) q[1];
sx q[1];
rz(-2.5320876) q[1];
sx q[1];
rz(-1.4896726) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4091858) q[0];
sx q[0];
rz(-0.86940765) q[0];
sx q[0];
rz(-2.5138445) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0428558) q[2];
sx q[2];
rz(-1.615881) q[2];
sx q[2];
rz(-2.2301607) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.5628964) q[1];
sx q[1];
rz(-1.7637196) q[1];
sx q[1];
rz(-1.3955411) q[1];
x q[2];
rz(1.0350111) q[3];
sx q[3];
rz(-1.1318558) q[3];
sx q[3];
rz(-0.5700232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5902517) q[2];
sx q[2];
rz(-2.6987023) q[2];
sx q[2];
rz(-1.7112973) q[2];
rz(-2.5643505) q[3];
sx q[3];
rz(-2.2539299) q[3];
sx q[3];
rz(-1.757471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5532613) q[0];
sx q[0];
rz(-1.7734779) q[0];
sx q[0];
rz(2.8531895) q[0];
rz(-0.53238955) q[1];
sx q[1];
rz(-0.45982292) q[1];
sx q[1];
rz(2.9945701) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0518236) q[0];
sx q[0];
rz(-0.71634403) q[0];
sx q[0];
rz(-2.1728974) q[0];
rz(-pi) q[1];
rz(0.27980079) q[2];
sx q[2];
rz(-2.8923312) q[2];
sx q[2];
rz(-1.3485497) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3210996) q[1];
sx q[1];
rz(-1.5309146) q[1];
sx q[1];
rz(0.95221968) q[1];
rz(0.57659984) q[3];
sx q[3];
rz(-1.7017662) q[3];
sx q[3];
rz(-0.031158202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0816575) q[2];
sx q[2];
rz(-2.7066878) q[2];
sx q[2];
rz(0.74404136) q[2];
rz(-2.384281) q[3];
sx q[3];
rz(-1.363874) q[3];
sx q[3];
rz(-1.7448759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1159146) q[0];
sx q[0];
rz(-1.0703351) q[0];
sx q[0];
rz(-1.0967789) q[0];
rz(-2.3241282) q[1];
sx q[1];
rz(-1.9348963) q[1];
sx q[1];
rz(2.5111326) q[1];
rz(-3.0030737) q[2];
sx q[2];
rz(-2.685683) q[2];
sx q[2];
rz(-1.4740623) q[2];
rz(0.13027262) q[3];
sx q[3];
rz(-2.1289005) q[3];
sx q[3];
rz(2.0990513) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
