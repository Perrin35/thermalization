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
rz(2.0107438) q[0];
sx q[0];
rz(-1.7717489) q[0];
sx q[0];
rz(-2.9983591) q[0];
rz(3.1205966) q[1];
sx q[1];
rz(3.723998) q[1];
sx q[1];
rz(5.8082685) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1256533) q[0];
sx q[0];
rz(-1.6467735) q[0];
sx q[0];
rz(-1.6267908) q[0];
rz(-pi) q[1];
rz(-2.5266746) q[2];
sx q[2];
rz(-1.5496829) q[2];
sx q[2];
rz(-0.17681618) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2432952) q[1];
sx q[1];
rz(-1.4615294) q[1];
sx q[1];
rz(-0.79078773) q[1];
rz(-pi) q[2];
rz(-1.1056312) q[3];
sx q[3];
rz(-2.0099927) q[3];
sx q[3];
rz(-0.10094563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1249866) q[2];
sx q[2];
rz(-1.6281444) q[2];
sx q[2];
rz(-2.5561257) q[2];
rz(-0.79750195) q[3];
sx q[3];
rz(-1.5471285) q[3];
sx q[3];
rz(2.7378979) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1095235) q[0];
sx q[0];
rz(-1.344498) q[0];
sx q[0];
rz(-2.4865785) q[0];
rz(0.0034927448) q[1];
sx q[1];
rz(-2.4066636) q[1];
sx q[1];
rz(0.014911501) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5302322) q[0];
sx q[0];
rz(-1.4900299) q[0];
sx q[0];
rz(1.5448891) q[0];
rz(-pi) q[1];
rz(-1.0833561) q[2];
sx q[2];
rz(-0.76448802) q[2];
sx q[2];
rz(-2.9202094) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7653489) q[1];
sx q[1];
rz(-2.2606008) q[1];
sx q[1];
rz(-2.6756555) q[1];
x q[2];
rz(-3.0767617) q[3];
sx q[3];
rz(-0.81063834) q[3];
sx q[3];
rz(-0.027934542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.94987503) q[2];
sx q[2];
rz(-1.3363375) q[2];
sx q[2];
rz(-2.8548262) q[2];
rz(3.1112572) q[3];
sx q[3];
rz(-1.560863) q[3];
sx q[3];
rz(-2.0386157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86998087) q[0];
sx q[0];
rz(-2.5219707) q[0];
sx q[0];
rz(1.9770589) q[0];
rz(-2.8097235) q[1];
sx q[1];
rz(-1.7121366) q[1];
sx q[1];
rz(-1.1480931) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5184421) q[0];
sx q[0];
rz(-1.352835) q[0];
sx q[0];
rz(0.017780546) q[0];
rz(-pi) q[1];
rz(2.0473061) q[2];
sx q[2];
rz(-1.4226573) q[2];
sx q[2];
rz(0.25853192) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.61473523) q[1];
sx q[1];
rz(-1.2678483) q[1];
sx q[1];
rz(-2.9399859) q[1];
rz(-pi) q[2];
rz(-1.8397736) q[3];
sx q[3];
rz(-0.1844953) q[3];
sx q[3];
rz(-2.6752145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1295192) q[2];
sx q[2];
rz(-1.4612863) q[2];
sx q[2];
rz(0.47895437) q[2];
rz(0.78914133) q[3];
sx q[3];
rz(-0.68500566) q[3];
sx q[3];
rz(1.3705378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-1.0399465) q[0];
sx q[0];
rz(-2.3532823) q[0];
sx q[0];
rz(-0.76860651) q[0];
rz(-2.6885314) q[1];
sx q[1];
rz(-2.1078347) q[1];
sx q[1];
rz(0.57910848) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.822389) q[0];
sx q[0];
rz(-1.3487848) q[0];
sx q[0];
rz(-0.25698203) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2007354) q[2];
sx q[2];
rz(-2.5343916) q[2];
sx q[2];
rz(3.0242357) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.107686) q[1];
sx q[1];
rz(-1.6158069) q[1];
sx q[1];
rz(2.0652524) q[1];
x q[2];
rz(-0.42570444) q[3];
sx q[3];
rz(-2.7085266) q[3];
sx q[3];
rz(-1.8830397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1733178) q[2];
sx q[2];
rz(-1.3710794) q[2];
sx q[2];
rz(-1.0260065) q[2];
rz(-0.46938986) q[3];
sx q[3];
rz(-1.0204851) q[3];
sx q[3];
rz(0.49218407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7475209) q[0];
sx q[0];
rz(-1.6265656) q[0];
sx q[0];
rz(-1.0006022) q[0];
rz(1.8383149) q[1];
sx q[1];
rz(-0.81221) q[1];
sx q[1];
rz(-2.2607048) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7422646) q[0];
sx q[0];
rz(-0.77631809) q[0];
sx q[0];
rz(-0.72350435) q[0];
rz(-2.4816031) q[2];
sx q[2];
rz(-1.3638328) q[2];
sx q[2];
rz(0.78850888) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.92915308) q[1];
sx q[1];
rz(-0.87599659) q[1];
sx q[1];
rz(1.9966182) q[1];
rz(-pi) q[2];
rz(0.88663574) q[3];
sx q[3];
rz(-0.58707844) q[3];
sx q[3];
rz(-2.1140703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.181043) q[2];
sx q[2];
rz(-1.1134104) q[2];
sx q[2];
rz(1.95365) q[2];
rz(-3.0175967) q[3];
sx q[3];
rz(-1.32722) q[3];
sx q[3];
rz(-2.7441062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7749629) q[0];
sx q[0];
rz(-0.22659817) q[0];
sx q[0];
rz(0.15283787) q[0];
rz(0.14640181) q[1];
sx q[1];
rz(-1.791879) q[1];
sx q[1];
rz(0.70405594) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25239736) q[0];
sx q[0];
rz(-1.5533698) q[0];
sx q[0];
rz(3.1118183) q[0];
rz(2.2302349) q[2];
sx q[2];
rz(-1.3040841) q[2];
sx q[2];
rz(-0.12202036) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4761281) q[1];
sx q[1];
rz(-1.4515467) q[1];
sx q[1];
rz(-0.075983451) q[1];
rz(-pi) q[2];
rz(1.1051768) q[3];
sx q[3];
rz(-2.4676968) q[3];
sx q[3];
rz(0.020581882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0600837) q[2];
sx q[2];
rz(-1.0241822) q[2];
sx q[2];
rz(-0.3788968) q[2];
rz(-1.0259519) q[3];
sx q[3];
rz(-0.70464269) q[3];
sx q[3];
rz(-2.0254693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41554552) q[0];
sx q[0];
rz(-1.3274095) q[0];
sx q[0];
rz(0.56336796) q[0];
rz(-0.4013966) q[1];
sx q[1];
rz(-1.0029663) q[1];
sx q[1];
rz(0.67515236) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9440985) q[0];
sx q[0];
rz(-1.5821274) q[0];
sx q[0];
rz(-1.5392153) q[0];
rz(-1.364351) q[2];
sx q[2];
rz(-2.1690302) q[2];
sx q[2];
rz(-0.75144671) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8842089) q[1];
sx q[1];
rz(-1.3335449) q[1];
sx q[1];
rz(-1.2409741) q[1];
rz(-pi) q[2];
rz(-2.3067857) q[3];
sx q[3];
rz(-1.2865598) q[3];
sx q[3];
rz(-2.4175197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1420765) q[2];
sx q[2];
rz(-2.2579305) q[2];
sx q[2];
rz(0.42259541) q[2];
rz(-2.8790867) q[3];
sx q[3];
rz(-1.0650977) q[3];
sx q[3];
rz(-2.0601823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-3.066684) q[0];
sx q[0];
rz(-0.96195641) q[0];
sx q[0];
rz(0.99895507) q[0];
rz(-1.9282438) q[1];
sx q[1];
rz(-1.1937001) q[1];
sx q[1];
rz(-2.8003069) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6741535) q[0];
sx q[0];
rz(-1.3786424) q[0];
sx q[0];
rz(-1.7056607) q[0];
x q[1];
rz(0.90088021) q[2];
sx q[2];
rz(-1.7770801) q[2];
sx q[2];
rz(-1.2717977) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5707376) q[1];
sx q[1];
rz(-2.8767817) q[1];
sx q[1];
rz(-0.45414208) q[1];
rz(-pi) q[2];
rz(-1.8496277) q[3];
sx q[3];
rz(-1.6270472) q[3];
sx q[3];
rz(-2.0771528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.21218941) q[2];
sx q[2];
rz(-1.5985649) q[2];
sx q[2];
rz(-2.8821778) q[2];
rz(2.8560396) q[3];
sx q[3];
rz(-0.67215896) q[3];
sx q[3];
rz(2.0937505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9015273) q[0];
sx q[0];
rz(-2.6492388) q[0];
sx q[0];
rz(2.0557851) q[0];
rz(-1.0049413) q[1];
sx q[1];
rz(-2.0045547) q[1];
sx q[1];
rz(-2.5850632) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96592007) q[0];
sx q[0];
rz(-1.5300599) q[0];
sx q[0];
rz(-1.7044742) q[0];
rz(0.6350216) q[2];
sx q[2];
rz(-0.94418923) q[2];
sx q[2];
rz(0.53900669) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.069935497) q[1];
sx q[1];
rz(-1.8462853) q[1];
sx q[1];
rz(-1.4749697) q[1];
rz(-pi) q[2];
rz(-1.0550761) q[3];
sx q[3];
rz(-2.0175126) q[3];
sx q[3];
rz(-1.8396256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.13653798) q[2];
sx q[2];
rz(-1.4056861) q[2];
sx q[2];
rz(1.8704937) q[2];
rz(-2.2942885) q[3];
sx q[3];
rz(-1.8151135) q[3];
sx q[3];
rz(-0.29720753) q[3];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6274881) q[0];
sx q[0];
rz(-1.4493554) q[0];
sx q[0];
rz(0.76747146) q[0];
rz(2.1766359) q[1];
sx q[1];
rz(-1.283353) q[1];
sx q[1];
rz(-2.8585785) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8022158) q[0];
sx q[0];
rz(-0.87516038) q[0];
sx q[0];
rz(2.8921769) q[0];
x q[1];
rz(-2.709952) q[2];
sx q[2];
rz(-1.5520763) q[2];
sx q[2];
rz(1.7016379) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9099847) q[1];
sx q[1];
rz(-0.63721133) q[1];
sx q[1];
rz(-0.63225327) q[1];
rz(-2.1238632) q[3];
sx q[3];
rz(-0.8700287) q[3];
sx q[3];
rz(-1.7297945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9877801) q[2];
sx q[2];
rz(-0.38975468) q[2];
sx q[2];
rz(-1.7392996) q[2];
rz(-0.66579372) q[3];
sx q[3];
rz(-1.8699162) q[3];
sx q[3];
rz(1.8388892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6233728) q[0];
sx q[0];
rz(-1.2517396) q[0];
sx q[0];
rz(-1.804833) q[0];
rz(1.5917336) q[1];
sx q[1];
rz(-1.5668329) q[1];
sx q[1];
rz(-0.11650539) q[1];
rz(0.62327173) q[2];
sx q[2];
rz(-2.4583271) q[2];
sx q[2];
rz(-1.5535977) q[2];
rz(0.79181087) q[3];
sx q[3];
rz(-1.4219501) q[3];
sx q[3];
rz(2.8590305) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
