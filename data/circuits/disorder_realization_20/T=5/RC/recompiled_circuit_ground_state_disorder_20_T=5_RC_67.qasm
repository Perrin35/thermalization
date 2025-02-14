OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5372758) q[0];
sx q[0];
rz(-0.24157) q[0];
sx q[0];
rz(0.33302745) q[0];
rz(-1.1355407) q[1];
sx q[1];
rz(3.9685213) q[1];
sx q[1];
rz(10.068738) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.052383311) q[0];
sx q[0];
rz(-2.174447) q[0];
sx q[0];
rz(-2.8103845) q[0];
rz(-pi) q[1];
rz(1.4139373) q[2];
sx q[2];
rz(-1.867395) q[2];
sx q[2];
rz(-0.083904249) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.79519479) q[1];
sx q[1];
rz(-1.9479706) q[1];
sx q[1];
rz(-1.8257797) q[1];
rz(0.34791874) q[3];
sx q[3];
rz(-1.9609465) q[3];
sx q[3];
rz(1.3044692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6551299) q[2];
sx q[2];
rz(-2.6772406) q[2];
sx q[2];
rz(0.74074024) q[2];
rz(-1.5517392) q[3];
sx q[3];
rz(-0.71151763) q[3];
sx q[3];
rz(-0.5927425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43316677) q[0];
sx q[0];
rz(-2.0080703) q[0];
sx q[0];
rz(-3.0246227) q[0];
rz(-0.50432694) q[1];
sx q[1];
rz(-1.3386936) q[1];
sx q[1];
rz(-0.28409827) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.741521) q[0];
sx q[0];
rz(-1.0062381) q[0];
sx q[0];
rz(-0.098640504) q[0];
rz(-pi) q[1];
rz(-2.1524736) q[2];
sx q[2];
rz(-2.2788413) q[2];
sx q[2];
rz(3.0281554) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5631905) q[1];
sx q[1];
rz(-3.0147073) q[1];
sx q[1];
rz(-1.2940426) q[1];
rz(1.8485214) q[3];
sx q[3];
rz(-1.3284995) q[3];
sx q[3];
rz(1.5577858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.32610193) q[2];
sx q[2];
rz(-0.88560605) q[2];
sx q[2];
rz(0.76078129) q[2];
rz(2.271999) q[3];
sx q[3];
rz(-1.2660916) q[3];
sx q[3];
rz(-2.6884955) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3608383) q[0];
sx q[0];
rz(-2.0212845) q[0];
sx q[0];
rz(-0.25217062) q[0];
rz(-1.699327) q[1];
sx q[1];
rz(-2.1863054) q[1];
sx q[1];
rz(2.7853277) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2748791) q[0];
sx q[0];
rz(-1.5769813) q[0];
sx q[0];
rz(1.6303819) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1596456) q[2];
sx q[2];
rz(-2.0817588) q[2];
sx q[2];
rz(1.1594349) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5327742) q[1];
sx q[1];
rz(-1.428481) q[1];
sx q[1];
rz(-2.2102093) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.0019508501) q[3];
sx q[3];
rz(-1.8271128) q[3];
sx q[3];
rz(1.0904097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0226655) q[2];
sx q[2];
rz(-1.7513195) q[2];
sx q[2];
rz(1.5125795) q[2];
rz(0.0096983612) q[3];
sx q[3];
rz(-1.2303979) q[3];
sx q[3];
rz(-2.5201918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6464722) q[0];
sx q[0];
rz(-2.5992114) q[0];
sx q[0];
rz(2.8908308) q[0];
rz(-1.3465025) q[1];
sx q[1];
rz(-0.52918068) q[1];
sx q[1];
rz(1.6983324) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.051606962) q[0];
sx q[0];
rz(-2.3653226) q[0];
sx q[0];
rz(1.2597741) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7378504) q[2];
sx q[2];
rz(-1.0761257) q[2];
sx q[2];
rz(0.97569377) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5574675) q[1];
sx q[1];
rz(-0.87601501) q[1];
sx q[1];
rz(-0.71895878) q[1];
rz(0.79899995) q[3];
sx q[3];
rz(-2.7741787) q[3];
sx q[3];
rz(1.1297117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.5556339) q[2];
sx q[2];
rz(-2.0363753) q[2];
sx q[2];
rz(-2.8982758) q[2];
rz(0.58468741) q[3];
sx q[3];
rz(-2.5644315) q[3];
sx q[3];
rz(3.1134591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1481767) q[0];
sx q[0];
rz(-2.7758444) q[0];
sx q[0];
rz(-2.962033) q[0];
rz(-1.2252294) q[1];
sx q[1];
rz(-1.7370677) q[1];
sx q[1];
rz(2.6833351) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4438547) q[0];
sx q[0];
rz(-2.8354044) q[0];
sx q[0];
rz(0.26060391) q[0];
x q[1];
rz(-2.8446537) q[2];
sx q[2];
rz(-2.3910948) q[2];
sx q[2];
rz(0.82210449) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1128487) q[1];
sx q[1];
rz(-0.15410559) q[1];
sx q[1];
rz(-1.0636368) q[1];
rz(-2.2716801) q[3];
sx q[3];
rz(-1.2906089) q[3];
sx q[3];
rz(0.82893574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7563584) q[2];
sx q[2];
rz(-2.7172654) q[2];
sx q[2];
rz(2.2198086) q[2];
rz(-1.2515757) q[3];
sx q[3];
rz(-1.4082963) q[3];
sx q[3];
rz(-3.0799227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.09403041) q[0];
sx q[0];
rz(-2.229409) q[0];
sx q[0];
rz(0.14866522) q[0];
rz(-2.8246763) q[1];
sx q[1];
rz(-0.47873679) q[1];
sx q[1];
rz(-2.5247578) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85108269) q[0];
sx q[0];
rz(-1.6491873) q[0];
sx q[0];
rz(0.061603994) q[0];
x q[1];
rz(-2.7814034) q[2];
sx q[2];
rz(-1.9860886) q[2];
sx q[2];
rz(-2.5915938) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3909797) q[1];
sx q[1];
rz(-1.3206078) q[1];
sx q[1];
rz(-0.98586086) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0265977) q[3];
sx q[3];
rz(-0.57145703) q[3];
sx q[3];
rz(0.28459099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2769015) q[2];
sx q[2];
rz(-1.428705) q[2];
sx q[2];
rz(-0.0083262715) q[2];
rz(-0.62638038) q[3];
sx q[3];
rz(-2.4255987) q[3];
sx q[3];
rz(3.1410419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.83547) q[0];
sx q[0];
rz(-1.1518814) q[0];
sx q[0];
rz(-0.48208153) q[0];
rz(0.029190633) q[1];
sx q[1];
rz(-1.8455285) q[1];
sx q[1];
rz(-0.76404244) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4493443) q[0];
sx q[0];
rz(-1.6381761) q[0];
sx q[0];
rz(1.9420366) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2305626) q[2];
sx q[2];
rz(-1.5507959) q[2];
sx q[2];
rz(2.3338712) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.58043874) q[1];
sx q[1];
rz(-2.517546) q[1];
sx q[1];
rz(1.8638205) q[1];
rz(-pi) q[2];
x q[2];
rz(1.589631) q[3];
sx q[3];
rz(-0.40098396) q[3];
sx q[3];
rz(0.31664059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.46334106) q[2];
sx q[2];
rz(-1.8222787) q[2];
sx q[2];
rz(-0.48416644) q[2];
rz(-0.68228996) q[3];
sx q[3];
rz(-2.4874918) q[3];
sx q[3];
rz(-0.13474034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.057137) q[0];
sx q[0];
rz(-2.3637922) q[0];
sx q[0];
rz(1.2917668) q[0];
rz(0.46503398) q[1];
sx q[1];
rz(-2.6227622) q[1];
sx q[1];
rz(2.8344287) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7503459) q[0];
sx q[0];
rz(-2.0758817) q[0];
sx q[0];
rz(-0.48663346) q[0];
x q[1];
rz(2.1910153) q[2];
sx q[2];
rz(-0.75846106) q[2];
sx q[2];
rz(1.4209227) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.41514709) q[1];
sx q[1];
rz(-1.4142805) q[1];
sx q[1];
rz(-0.99516258) q[1];
rz(-0.62509663) q[3];
sx q[3];
rz(-1.5872883) q[3];
sx q[3];
rz(-2.0979289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9610567) q[2];
sx q[2];
rz(-2.7013216) q[2];
sx q[2];
rz(-0.72009909) q[2];
rz(-0.85047203) q[3];
sx q[3];
rz(-1.1668147) q[3];
sx q[3];
rz(-1.7299962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20486031) q[0];
sx q[0];
rz(-2.9731049) q[0];
sx q[0];
rz(0.70814651) q[0];
rz(-1.5180961) q[1];
sx q[1];
rz(-2.203439) q[1];
sx q[1];
rz(0.35619563) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30763876) q[0];
sx q[0];
rz(-2.4705268) q[0];
sx q[0];
rz(-2.4249581) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6004066) q[2];
sx q[2];
rz(-2.1884754) q[2];
sx q[2];
rz(0.35770616) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.5070008) q[1];
sx q[1];
rz(-2.2193877) q[1];
sx q[1];
rz(-1.6084987) q[1];
x q[2];
rz(-2.735504) q[3];
sx q[3];
rz(-2.0022939) q[3];
sx q[3];
rz(0.0270947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0922962) q[2];
sx q[2];
rz(-2.3342817) q[2];
sx q[2];
rz(-0.88579196) q[2];
rz(-0.93585912) q[3];
sx q[3];
rz(-1.9610619) q[3];
sx q[3];
rz(0.22127557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43779272) q[0];
sx q[0];
rz(-3.0942823) q[0];
sx q[0];
rz(1.6920775) q[0];
rz(2.119078) q[1];
sx q[1];
rz(-0.48373628) q[1];
sx q[1];
rz(-1.6922916) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81240772) q[0];
sx q[0];
rz(-1.9827166) q[0];
sx q[0];
rz(-3.0874041) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76646508) q[2];
sx q[2];
rz(-2.577707) q[2];
sx q[2];
rz(-1.5001129) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.76043452) q[1];
sx q[1];
rz(-1.0236003) q[1];
sx q[1];
rz(-0.24914279) q[1];
x q[2];
rz(2.0725771) q[3];
sx q[3];
rz(-0.3007362) q[3];
sx q[3];
rz(-1.6274483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.25541043) q[2];
sx q[2];
rz(-1.9516727) q[2];
sx q[2];
rz(2.4667242) q[2];
rz(0.97149649) q[3];
sx q[3];
rz(-0.70236218) q[3];
sx q[3];
rz(1.6500047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3424727) q[0];
sx q[0];
rz(-1.5958888) q[0];
sx q[0];
rz(2.2842443) q[0];
rz(0.2324066) q[1];
sx q[1];
rz(-2.1288165) q[1];
sx q[1];
rz(-1.7850599) q[1];
rz(-0.92844456) q[2];
sx q[2];
rz(-2.5088163) q[2];
sx q[2];
rz(-2.8586254) q[2];
rz(-1.6583937) q[3];
sx q[3];
rz(-2.4110553) q[3];
sx q[3];
rz(0.72881107) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
