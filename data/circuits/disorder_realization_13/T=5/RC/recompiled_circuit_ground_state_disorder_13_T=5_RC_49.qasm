OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.00081113022) q[0];
sx q[0];
rz(-0.82511628) q[0];
sx q[0];
rz(-2.7197279) q[0];
rz(2.3946664) q[1];
sx q[1];
rz(-0.59023017) q[1];
sx q[1];
rz(-2.3205369) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79347389) q[0];
sx q[0];
rz(-3.0464026) q[0];
sx q[0];
rz(-0.4918672) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9157588) q[2];
sx q[2];
rz(-1.7712812) q[2];
sx q[2];
rz(-0.12910138) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3154308) q[1];
sx q[1];
rz(-1.9524442) q[1];
sx q[1];
rz(2.8355469) q[1];
rz(-1.2263359) q[3];
sx q[3];
rz(-1.4446961) q[3];
sx q[3];
rz(0.9094224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7464298) q[2];
sx q[2];
rz(-0.91150993) q[2];
sx q[2];
rz(0.54473031) q[2];
rz(1.4970477) q[3];
sx q[3];
rz(-1.112273) q[3];
sx q[3];
rz(-1.8018319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4561975) q[0];
sx q[0];
rz(-2.39769) q[0];
sx q[0];
rz(0.5734545) q[0];
rz(-0.042512976) q[1];
sx q[1];
rz(-1.7492234) q[1];
sx q[1];
rz(1.858985) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61717466) q[0];
sx q[0];
rz(-1.5502009) q[0];
sx q[0];
rz(-3.0153946) q[0];
x q[1];
rz(1.3147215) q[2];
sx q[2];
rz(-1.2856069) q[2];
sx q[2];
rz(1.3392753) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.32795045) q[1];
sx q[1];
rz(-0.54681603) q[1];
sx q[1];
rz(-0.27346404) q[1];
x q[2];
rz(2.1759791) q[3];
sx q[3];
rz(-2.11016) q[3];
sx q[3];
rz(0.920528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.87806988) q[2];
sx q[2];
rz(-1.9922549) q[2];
sx q[2];
rz(2.6980706) q[2];
rz(1.5524607) q[3];
sx q[3];
rz(-2.7510721) q[3];
sx q[3];
rz(-0.21903567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1294915) q[0];
sx q[0];
rz(-0.87915593) q[0];
sx q[0];
rz(-2.7398859) q[0];
rz(1.4178287) q[1];
sx q[1];
rz(-2.6938853) q[1];
sx q[1];
rz(-2.0254501) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6594769) q[0];
sx q[0];
rz(-1.2667313) q[0];
sx q[0];
rz(-0.67993865) q[0];
rz(-pi) q[1];
rz(0.1537519) q[2];
sx q[2];
rz(-2.2567686) q[2];
sx q[2];
rz(-1.5123715) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1895208) q[1];
sx q[1];
rz(-2.4671989) q[1];
sx q[1];
rz(2.2728069) q[1];
x q[2];
rz(0.54858923) q[3];
sx q[3];
rz(-1.5058796) q[3];
sx q[3];
rz(1.3260076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9125354) q[2];
sx q[2];
rz(-2.8450862) q[2];
sx q[2];
rz(-2.1336446) q[2];
rz(-1.6352765) q[3];
sx q[3];
rz(-1.1791891) q[3];
sx q[3];
rz(-0.87872163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2849543) q[0];
sx q[0];
rz(-3.1143739) q[0];
sx q[0];
rz(0.83754367) q[0];
rz(-1.2614999) q[1];
sx q[1];
rz(-1.7439758) q[1];
sx q[1];
rz(-1.8998442) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57056242) q[0];
sx q[0];
rz(-1.116291) q[0];
sx q[0];
rz(-1.8733371) q[0];
rz(-pi) q[1];
rz(2.9292232) q[2];
sx q[2];
rz(-0.97281934) q[2];
sx q[2];
rz(0.47456196) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.82363331) q[1];
sx q[1];
rz(-1.1640062) q[1];
sx q[1];
rz(0.49655668) q[1];
rz(-pi) q[2];
rz(-0.37184985) q[3];
sx q[3];
rz(-2.1193373) q[3];
sx q[3];
rz(-1.8601102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1602829) q[2];
sx q[2];
rz(-2.1872988) q[2];
sx q[2];
rz(2.1295638) q[2];
rz(-0.99803287) q[3];
sx q[3];
rz(-2.4055552) q[3];
sx q[3];
rz(1.6691104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0618458) q[0];
sx q[0];
rz(-2.340402) q[0];
sx q[0];
rz(2.8629942) q[0];
rz(-2.8246236) q[1];
sx q[1];
rz(-1.3895037) q[1];
sx q[1];
rz(-0.46151361) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7636895) q[0];
sx q[0];
rz(-0.17556854) q[0];
sx q[0];
rz(-1.8493091) q[0];
rz(-pi) q[1];
x q[1];
rz(0.95038173) q[2];
sx q[2];
rz(-1.912552) q[2];
sx q[2];
rz(2.1555962) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0089142) q[1];
sx q[1];
rz(-1.8092833) q[1];
sx q[1];
rz(-0.096367457) q[1];
rz(3.0449681) q[3];
sx q[3];
rz(-1.0276252) q[3];
sx q[3];
rz(-2.4718049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3658112) q[2];
sx q[2];
rz(-1.3373988) q[2];
sx q[2];
rz(-2.865045) q[2];
rz(0.60025275) q[3];
sx q[3];
rz(-2.4252031) q[3];
sx q[3];
rz(-2.6071809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9944331) q[0];
sx q[0];
rz(-0.58013791) q[0];
sx q[0];
rz(-0.68122) q[0];
rz(2.6874806) q[1];
sx q[1];
rz(-1.157016) q[1];
sx q[1];
rz(-2.5195794) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3709741) q[0];
sx q[0];
rz(-1.7683836) q[0];
sx q[0];
rz(1.4739735) q[0];
rz(-pi) q[1];
x q[1];
rz(0.95152835) q[2];
sx q[2];
rz(-1.027809) q[2];
sx q[2];
rz(-1.1752626) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6890735) q[1];
sx q[1];
rz(-1.9770685) q[1];
sx q[1];
rz(-2.3157332) q[1];
rz(-pi) q[2];
rz(-2.8735749) q[3];
sx q[3];
rz(-1.9081731) q[3];
sx q[3];
rz(1.5609139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6260208) q[2];
sx q[2];
rz(-0.61352789) q[2];
sx q[2];
rz(-2.386509) q[2];
rz(-3.0355022) q[3];
sx q[3];
rz(-1.875149) q[3];
sx q[3];
rz(0.46060002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.047091529) q[0];
sx q[0];
rz(-0.61151183) q[0];
sx q[0];
rz(3.0754572) q[0];
rz(0.051008929) q[1];
sx q[1];
rz(-2.3936733) q[1];
sx q[1];
rz(1.8775108) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8569778) q[0];
sx q[0];
rz(-0.98941411) q[0];
sx q[0];
rz(-0.435243) q[0];
rz(-pi) q[1];
rz(-0.35392742) q[2];
sx q[2];
rz(-2.3945237) q[2];
sx q[2];
rz(-0.22176556) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7552774) q[1];
sx q[1];
rz(-1.9615972) q[1];
sx q[1];
rz(-2.3936207) q[1];
rz(-pi) q[2];
rz(-1.966428) q[3];
sx q[3];
rz(-1.9744748) q[3];
sx q[3];
rz(0.56433557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.85822004) q[2];
sx q[2];
rz(-1.1274575) q[2];
sx q[2];
rz(2.9952725) q[2];
rz(-2.537651) q[3];
sx q[3];
rz(-0.54655176) q[3];
sx q[3];
rz(-1.4124136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42501763) q[0];
sx q[0];
rz(-1.1973493) q[0];
sx q[0];
rz(0.091751598) q[0];
rz(-3.0692406) q[1];
sx q[1];
rz(-1.3183343) q[1];
sx q[1];
rz(2.4718557) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90988542) q[0];
sx q[0];
rz(-0.32389952) q[0];
sx q[0];
rz(-0.43816752) q[0];
rz(-pi) q[1];
rz(0.40939949) q[2];
sx q[2];
rz(-1.4509307) q[2];
sx q[2];
rz(1.9425962) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.200676) q[1];
sx q[1];
rz(-1.426061) q[1];
sx q[1];
rz(2.8334068) q[1];
rz(-pi) q[2];
rz(-2.9961139) q[3];
sx q[3];
rz(-1.1324185) q[3];
sx q[3];
rz(-2.9117025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3503795) q[2];
sx q[2];
rz(-0.71724856) q[2];
sx q[2];
rz(2.4199602) q[2];
rz(-0.77389884) q[3];
sx q[3];
rz(-0.98267233) q[3];
sx q[3];
rz(2.6017046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4197107) q[0];
sx q[0];
rz(-1.1996491) q[0];
sx q[0];
rz(-2.6019959) q[0];
rz(-2.3609912) q[1];
sx q[1];
rz(-1.8949948) q[1];
sx q[1];
rz(-2.7194729) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0142067) q[0];
sx q[0];
rz(-2.2932568) q[0];
sx q[0];
rz(2.3425472) q[0];
rz(-3.0303553) q[2];
sx q[2];
rz(-2.6468177) q[2];
sx q[2];
rz(0.24559427) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5663492) q[1];
sx q[1];
rz(-0.48108992) q[1];
sx q[1];
rz(-1.653137) q[1];
rz(-pi) q[2];
rz(-0.5886505) q[3];
sx q[3];
rz(-1.5331309) q[3];
sx q[3];
rz(-1.879694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0427986) q[2];
sx q[2];
rz(-2.1197987) q[2];
sx q[2];
rz(-0.54879028) q[2];
rz(2.9889066) q[3];
sx q[3];
rz(-2.1137674) q[3];
sx q[3];
rz(-0.71167439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3720836) q[0];
sx q[0];
rz(-2.7476269) q[0];
sx q[0];
rz(2.5373996) q[0];
rz(1.4671885) q[1];
sx q[1];
rz(-1.3507651) q[1];
sx q[1];
rz(-1.9942572) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7580099) q[0];
sx q[0];
rz(-2.1910163) q[0];
sx q[0];
rz(2.0877286) q[0];
rz(-2.7962923) q[2];
sx q[2];
rz(-0.7658813) q[2];
sx q[2];
rz(1.6102546) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1318329) q[1];
sx q[1];
rz(-1.5861675) q[1];
sx q[1];
rz(-2.6587376) q[1];
rz(-pi) q[2];
rz(1.9705521) q[3];
sx q[3];
rz(-1.2467524) q[3];
sx q[3];
rz(-0.43982616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4597822) q[2];
sx q[2];
rz(-1.1638389) q[2];
sx q[2];
rz(-0.60010827) q[2];
rz(2.4188304) q[3];
sx q[3];
rz(-0.99812752) q[3];
sx q[3];
rz(2.7616937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73846524) q[0];
sx q[0];
rz(-1.6275788) q[0];
sx q[0];
rz(1.7049261) q[0];
rz(-1.5858831) q[1];
sx q[1];
rz(-1.7441505) q[1];
sx q[1];
rz(2.2081262) q[1];
rz(-2.2384833) q[2];
sx q[2];
rz(-1.9676859) q[2];
sx q[2];
rz(-3.0953593) q[2];
rz(-2.4248471) q[3];
sx q[3];
rz(-2.8195753) q[3];
sx q[3];
rz(-1.4618518) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
