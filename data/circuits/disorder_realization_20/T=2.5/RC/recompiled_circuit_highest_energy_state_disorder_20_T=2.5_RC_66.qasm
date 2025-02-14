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
rz(2.4169154) q[0];
sx q[0];
rz(-1.3114572) q[0];
sx q[0];
rz(-2.3367982) q[0];
rz(0.59504396) q[1];
sx q[1];
rz(-2.9437328) q[1];
sx q[1];
rz(-1.9207538) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9210584) q[0];
sx q[0];
rz(-3.0166114) q[0];
sx q[0];
rz(0.28470953) q[0];
x q[1];
rz(1.7096504) q[2];
sx q[2];
rz(-2.029691) q[2];
sx q[2];
rz(0.51990055) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8960305) q[1];
sx q[1];
rz(-2.3136825) q[1];
sx q[1];
rz(-1.6309392) q[1];
x q[2];
rz(1.4335459) q[3];
sx q[3];
rz(-0.94232925) q[3];
sx q[3];
rz(-2.8253205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7120984) q[2];
sx q[2];
rz(-0.50769371) q[2];
sx q[2];
rz(1.8053767) q[2];
rz(1.3974238) q[3];
sx q[3];
rz(-1.5066527) q[3];
sx q[3];
rz(1.5857182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76032388) q[0];
sx q[0];
rz(-1.5070494) q[0];
sx q[0];
rz(2.4632578) q[0];
rz(0.61596576) q[1];
sx q[1];
rz(-2.1149642) q[1];
sx q[1];
rz(-0.14332992) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5943439) q[0];
sx q[0];
rz(-1.4137725) q[0];
sx q[0];
rz(-1.2444522) q[0];
x q[1];
rz(0.11516659) q[2];
sx q[2];
rz(-0.85571849) q[2];
sx q[2];
rz(0.63181704) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.60494938) q[1];
sx q[1];
rz(-2.6338913) q[1];
sx q[1];
rz(-0.097670876) q[1];
rz(-2.7160104) q[3];
sx q[3];
rz(-1.2493361) q[3];
sx q[3];
rz(1.5150573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4914638) q[2];
sx q[2];
rz(-2.3626707) q[2];
sx q[2];
rz(1.9219363) q[2];
rz(-2.8236112) q[3];
sx q[3];
rz(-2.3173083) q[3];
sx q[3];
rz(2.4021751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3984482) q[0];
sx q[0];
rz(-0.92107451) q[0];
sx q[0];
rz(0.63968023) q[0];
rz(-1.8094874) q[1];
sx q[1];
rz(-0.94146856) q[1];
sx q[1];
rz(-2.926362) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9063961) q[0];
sx q[0];
rz(-2.0590889) q[0];
sx q[0];
rz(-1.8129691) q[0];
x q[1];
rz(1.8551183) q[2];
sx q[2];
rz(-1.5657565) q[2];
sx q[2];
rz(-1.3748295) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.58956026) q[1];
sx q[1];
rz(-1.5812687) q[1];
sx q[1];
rz(-1.9661436) q[1];
x q[2];
rz(-2.027719) q[3];
sx q[3];
rz(-1.6789241) q[3];
sx q[3];
rz(2.7939452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8362063) q[2];
sx q[2];
rz(-2.6658194) q[2];
sx q[2];
rz(2.9212908) q[2];
rz(0.78271714) q[3];
sx q[3];
rz(-1.3681151) q[3];
sx q[3];
rz(2.9960184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.7102605) q[0];
sx q[0];
rz(-2.1556957) q[0];
sx q[0];
rz(1.8248935) q[0];
rz(0.45007625) q[1];
sx q[1];
rz(-1.7066259) q[1];
sx q[1];
rz(-3.0302474) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25791439) q[0];
sx q[0];
rz(-2.2131569) q[0];
sx q[0];
rz(-1.2538373) q[0];
rz(-pi) q[1];
rz(3.0527332) q[2];
sx q[2];
rz(-1.2496762) q[2];
sx q[2];
rz(-0.027983657) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7584472) q[1];
sx q[1];
rz(-1.0186334) q[1];
sx q[1];
rz(-2.8124185) q[1];
x q[2];
rz(-2.4347859) q[3];
sx q[3];
rz(-0.68663678) q[3];
sx q[3];
rz(-1.6345018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.7340118) q[2];
sx q[2];
rz(-2.0110726) q[2];
sx q[2];
rz(-0.13047516) q[2];
rz(-2.981251) q[3];
sx q[3];
rz(-0.66157833) q[3];
sx q[3];
rz(-2.9714238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.570356) q[0];
sx q[0];
rz(-0.2061051) q[0];
sx q[0];
rz(-0.11827949) q[0];
rz(2.2453399) q[1];
sx q[1];
rz(-1.4786485) q[1];
sx q[1];
rz(0.55312696) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.129707) q[0];
sx q[0];
rz(-1.0965523) q[0];
sx q[0];
rz(-1.2007196) q[0];
rz(-2.4055355) q[2];
sx q[2];
rz(-1.9070574) q[2];
sx q[2];
rz(-2.2081809) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.70870465) q[1];
sx q[1];
rz(-2.1498393) q[1];
sx q[1];
rz(2.9396073) q[1];
x q[2];
rz(1.0712264) q[3];
sx q[3];
rz(-2.4132846) q[3];
sx q[3];
rz(2.0508931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7913738) q[2];
sx q[2];
rz(-0.72177902) q[2];
sx q[2];
rz(2.1283894) q[2];
rz(0.31747097) q[3];
sx q[3];
rz(-1.443913) q[3];
sx q[3];
rz(-2.1875994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1182227) q[0];
sx q[0];
rz(-2.2582) q[0];
sx q[0];
rz(-0.50514847) q[0];
rz(-2.743424) q[1];
sx q[1];
rz(-1.265637) q[1];
sx q[1];
rz(1.8123951) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1627462) q[0];
sx q[0];
rz(-1.1067451) q[0];
sx q[0];
rz(-0.18382182) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2114491) q[2];
sx q[2];
rz(-2.076175) q[2];
sx q[2];
rz(-2.0063248) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8230629) q[1];
sx q[1];
rz(-1.6803871) q[1];
sx q[1];
rz(-1.0062045) q[1];
rz(0.81327) q[3];
sx q[3];
rz(-2.7754042) q[3];
sx q[3];
rz(2.1048819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.087611) q[2];
sx q[2];
rz(-1.8113281) q[2];
sx q[2];
rz(2.3114253) q[2];
rz(1.6603445) q[3];
sx q[3];
rz(-2.0989428) q[3];
sx q[3];
rz(1.164485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2723349) q[0];
sx q[0];
rz(-2.2255958) q[0];
sx q[0];
rz(1.9328312) q[0];
rz(0.2441949) q[1];
sx q[1];
rz(-1.1714275) q[1];
sx q[1];
rz(0.29921439) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1580762) q[0];
sx q[0];
rz(-1.354748) q[0];
sx q[0];
rz(-3.1147577) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60986477) q[2];
sx q[2];
rz(-2.0195896) q[2];
sx q[2];
rz(3.0013709) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0639887) q[1];
sx q[1];
rz(-2.4875459) q[1];
sx q[1];
rz(1.4908423) q[1];
rz(-pi) q[2];
rz(0.47513606) q[3];
sx q[3];
rz(-2.3804602) q[3];
sx q[3];
rz(-0.5055529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9485335) q[2];
sx q[2];
rz(-1.7391354) q[2];
sx q[2];
rz(-2.194727) q[2];
rz(-1.7638505) q[3];
sx q[3];
rz(-0.54072127) q[3];
sx q[3];
rz(-0.64259678) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68719012) q[0];
sx q[0];
rz(-0.43455046) q[0];
sx q[0];
rz(2.6026671) q[0];
rz(2.9680805) q[1];
sx q[1];
rz(-0.75650802) q[1];
sx q[1];
rz(0.38421806) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.980775) q[0];
sx q[0];
rz(-1.4133102) q[0];
sx q[0];
rz(-2.6927267) q[0];
rz(0.27663934) q[2];
sx q[2];
rz(-0.62046548) q[2];
sx q[2];
rz(2.7850604) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9335971) q[1];
sx q[1];
rz(-1.5969445) q[1];
sx q[1];
rz(1.5118062) q[1];
rz(-pi) q[2];
rz(1.5269856) q[3];
sx q[3];
rz(-2.5871944) q[3];
sx q[3];
rz(1.2498433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3197799) q[2];
sx q[2];
rz(-0.76924789) q[2];
sx q[2];
rz(-0.5963076) q[2];
rz(-1.0373621) q[3];
sx q[3];
rz(-0.77330247) q[3];
sx q[3];
rz(-1.1843225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3768815) q[0];
sx q[0];
rz(-0.75279555) q[0];
sx q[0];
rz(2.9994614) q[0];
rz(1.1171974) q[1];
sx q[1];
rz(-1.6551207) q[1];
sx q[1];
rz(-1.2275009) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9949029) q[0];
sx q[0];
rz(-1.0953383) q[0];
sx q[0];
rz(1.2593996) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4551444) q[2];
sx q[2];
rz(-0.1384379) q[2];
sx q[2];
rz(-1.1756681) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.092526768) q[1];
sx q[1];
rz(-2.5124308) q[1];
sx q[1];
rz(-1.4607304) q[1];
rz(0.82955011) q[3];
sx q[3];
rz(-1.2365474) q[3];
sx q[3];
rz(3.0544508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.050798945) q[2];
sx q[2];
rz(-2.8992081) q[2];
sx q[2];
rz(0.79831115) q[2];
rz(1.5797041) q[3];
sx q[3];
rz(-1.3825994) q[3];
sx q[3];
rz(-0.17670512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0853737) q[0];
sx q[0];
rz(-2.4985785) q[0];
sx q[0];
rz(-3.1363078) q[0];
rz(0.082848631) q[1];
sx q[1];
rz(-2.0382035) q[1];
sx q[1];
rz(2.8740035) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8092958) q[0];
sx q[0];
rz(-2.3837261) q[0];
sx q[0];
rz(2.6790455) q[0];
rz(-1.9584974) q[2];
sx q[2];
rz(-1.299721) q[2];
sx q[2];
rz(0.77927206) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6050075) q[1];
sx q[1];
rz(-1.3571902) q[1];
sx q[1];
rz(-1.8940635) q[1];
rz(-0.4739828) q[3];
sx q[3];
rz(-1.1055531) q[3];
sx q[3];
rz(1.3657692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.573632) q[2];
sx q[2];
rz(-2.6601514) q[2];
sx q[2];
rz(1.2130223) q[2];
rz(1.8968449) q[3];
sx q[3];
rz(-0.48143482) q[3];
sx q[3];
rz(-1.1904233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.5382814) q[0];
sx q[0];
rz(-0.9700226) q[0];
sx q[0];
rz(0.097401311) q[0];
rz(3.1311323) q[1];
sx q[1];
rz(-0.86581007) q[1];
sx q[1];
rz(2.5444358) q[1];
rz(0.14870208) q[2];
sx q[2];
rz(-0.54051334) q[2];
sx q[2];
rz(1.9220026) q[2];
rz(-0.17531323) q[3];
sx q[3];
rz(-0.98876035) q[3];
sx q[3];
rz(0.50553346) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
