OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4133889) q[0];
sx q[0];
rz(-1.1336741) q[0];
sx q[0];
rz(1.5925621) q[0];
rz(1.6917317) q[1];
sx q[1];
rz(5.6258968) q[1];
sx q[1];
rz(13.110553) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9323174) q[0];
sx q[0];
rz(-1.6074751) q[0];
sx q[0];
rz(3.0741312) q[0];
rz(-3.0402096) q[2];
sx q[2];
rz(-1.6147436) q[2];
sx q[2];
rz(-1.4506538) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.45317263) q[1];
sx q[1];
rz(-1.9721713) q[1];
sx q[1];
rz(-2.32262) q[1];
rz(-pi) q[2];
rz(1.5531494) q[3];
sx q[3];
rz(-1.5994161) q[3];
sx q[3];
rz(-2.0538405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0698174) q[2];
sx q[2];
rz(-1.2640307) q[2];
sx q[2];
rz(-1.7791746) q[2];
rz(-3.1133364) q[3];
sx q[3];
rz(-1.3794206) q[3];
sx q[3];
rz(-0.79022592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4966999) q[0];
sx q[0];
rz(-2.4364478) q[0];
sx q[0];
rz(-1.9702966) q[0];
rz(2.9303739) q[1];
sx q[1];
rz(-0.44208458) q[1];
sx q[1];
rz(-1.404095) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068721213) q[0];
sx q[0];
rz(-2.255548) q[0];
sx q[0];
rz(2.1268334) q[0];
rz(-1.7369466) q[2];
sx q[2];
rz(-2.0478021) q[2];
sx q[2];
rz(-2.0376861) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.029570015) q[1];
sx q[1];
rz(-0.28407541) q[1];
sx q[1];
rz(3.0594538) q[1];
x q[2];
rz(-1.4349798) q[3];
sx q[3];
rz(-1.0944301) q[3];
sx q[3];
rz(-1.3729031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.30963787) q[2];
sx q[2];
rz(-1.4761304) q[2];
sx q[2];
rz(-0.94397604) q[2];
rz(0.55654636) q[3];
sx q[3];
rz(-0.49749938) q[3];
sx q[3];
rz(-1.336162) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29207644) q[0];
sx q[0];
rz(-2.3777666) q[0];
sx q[0];
rz(-1.7180432) q[0];
rz(0.81958333) q[1];
sx q[1];
rz(-2.3312566) q[1];
sx q[1];
rz(2.5779285) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.992313) q[0];
sx q[0];
rz(-2.5631319) q[0];
sx q[0];
rz(-3.0119386) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5416449) q[2];
sx q[2];
rz(-1.3834582) q[2];
sx q[2];
rz(1.3647321) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0020395) q[1];
sx q[1];
rz(-0.66322749) q[1];
sx q[1];
rz(-2.1173649) q[1];
rz(-pi) q[2];
rz(-1.9824969) q[3];
sx q[3];
rz(-0.99038306) q[3];
sx q[3];
rz(1.2702277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6391969) q[2];
sx q[2];
rz(-1.1221308) q[2];
sx q[2];
rz(-1.0602661) q[2];
rz(-3.0660196) q[3];
sx q[3];
rz(-0.85791701) q[3];
sx q[3];
rz(0.11463541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8183427) q[0];
sx q[0];
rz(-0.30775726) q[0];
sx q[0];
rz(-0.19317214) q[0];
rz(-1.5441783) q[1];
sx q[1];
rz(-1.1491821) q[1];
sx q[1];
rz(-0.65778041) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5415216) q[0];
sx q[0];
rz(-0.17734781) q[0];
sx q[0];
rz(-2.970201) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1409608) q[2];
sx q[2];
rz(-1.2865215) q[2];
sx q[2];
rz(2.5179799) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.66637908) q[1];
sx q[1];
rz(-0.5961601) q[1];
sx q[1];
rz(-2.0344884) q[1];
x q[2];
rz(0.25810453) q[3];
sx q[3];
rz(-1.8542284) q[3];
sx q[3];
rz(-2.5131445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0956991) q[2];
sx q[2];
rz(-0.4102439) q[2];
sx q[2];
rz(2.0945385) q[2];
rz(2.0783157) q[3];
sx q[3];
rz(-1.1626817) q[3];
sx q[3];
rz(-1.261238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3445774) q[0];
sx q[0];
rz(-2.8224967) q[0];
sx q[0];
rz(1.0239333) q[0];
rz(-0.78760415) q[1];
sx q[1];
rz(-1.5658295) q[1];
sx q[1];
rz(2.5147298) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1889362) q[0];
sx q[0];
rz(-0.034010012) q[0];
sx q[0];
rz(1.0556428) q[0];
rz(-2.3847694) q[2];
sx q[2];
rz(-2.7177817) q[2];
sx q[2];
rz(-0.28048453) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6080731) q[1];
sx q[1];
rz(-2.7058209) q[1];
sx q[1];
rz(-0.79969745) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8623779) q[3];
sx q[3];
rz(-0.40816187) q[3];
sx q[3];
rz(0.51418958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.91288599) q[2];
sx q[2];
rz(-1.9084946) q[2];
sx q[2];
rz(-2.1234925) q[2];
rz(-0.61156887) q[3];
sx q[3];
rz(-1.3221778) q[3];
sx q[3];
rz(2.1900246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
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
rz(1.4488895) q[0];
sx q[0];
rz(-1.792181) q[0];
sx q[0];
rz(-1.942379) q[0];
rz(-1.9723643) q[1];
sx q[1];
rz(-2.0904082) q[1];
sx q[1];
rz(0.32454023) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6231461) q[0];
sx q[0];
rz(-2.0707957) q[0];
sx q[0];
rz(-0.48796939) q[0];
x q[1];
rz(2.4915699) q[2];
sx q[2];
rz(-2.2215002) q[2];
sx q[2];
rz(2.0536154) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2734087) q[1];
sx q[1];
rz(-1.5460099) q[1];
sx q[1];
rz(2.1048057) q[1];
x q[2];
rz(-1.2732029) q[3];
sx q[3];
rz(-1.9801567) q[3];
sx q[3];
rz(-2.5149432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4679608) q[2];
sx q[2];
rz(-1.8362074) q[2];
sx q[2];
rz(1.8661873) q[2];
rz(-1.0547137) q[3];
sx q[3];
rz(-0.79189363) q[3];
sx q[3];
rz(1.595343) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4642898) q[0];
sx q[0];
rz(-1.807656) q[0];
sx q[0];
rz(1.4087079) q[0];
rz(-2.6858221) q[1];
sx q[1];
rz(-2.9401638) q[1];
sx q[1];
rz(1.2021525) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5613272) q[0];
sx q[0];
rz(-0.94764493) q[0];
sx q[0];
rz(-0.49669637) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1832778) q[2];
sx q[2];
rz(-2.2829208) q[2];
sx q[2];
rz(0.14856635) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.251293) q[1];
sx q[1];
rz(-1.4874465) q[1];
sx q[1];
rz(-3.0492196) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2644516) q[3];
sx q[3];
rz(-1.3925941) q[3];
sx q[3];
rz(1.7314535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1414286) q[2];
sx q[2];
rz(-1.2951853) q[2];
sx q[2];
rz(-2.2953575) q[2];
rz(-0.49514654) q[3];
sx q[3];
rz(-2.4954002) q[3];
sx q[3];
rz(1.5406066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071035944) q[0];
sx q[0];
rz(-1.8510011) q[0];
sx q[0];
rz(2.0284247) q[0];
rz(2.44599) q[1];
sx q[1];
rz(-0.38989392) q[1];
sx q[1];
rz(1.6092469) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42513645) q[0];
sx q[0];
rz(-1.3279337) q[0];
sx q[0];
rz(-2.6572685) q[0];
rz(-2.3959827) q[2];
sx q[2];
rz(-0.62586212) q[2];
sx q[2];
rz(1.9288043) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.53837) q[1];
sx q[1];
rz(-2.6909628) q[1];
sx q[1];
rz(-1.5249114) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.60204695) q[3];
sx q[3];
rz(-1.6685408) q[3];
sx q[3];
rz(-2.5412113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9939076) q[2];
sx q[2];
rz(-2.9116178) q[2];
sx q[2];
rz(0.85419401) q[2];
rz(-1.9130075) q[3];
sx q[3];
rz(-1.5649786) q[3];
sx q[3];
rz(1.6555697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5937186) q[0];
sx q[0];
rz(-2.6429208) q[0];
sx q[0];
rz(-0.6643995) q[0];
rz(1.9539072) q[1];
sx q[1];
rz(-0.70078754) q[1];
sx q[1];
rz(-1.2845576) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4391543) q[0];
sx q[0];
rz(-0.87411532) q[0];
sx q[0];
rz(-0.37581635) q[0];
x q[1];
rz(-0.2662439) q[2];
sx q[2];
rz(-1.9387445) q[2];
sx q[2];
rz(-3.0014696) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.38332332) q[1];
sx q[1];
rz(-0.82793923) q[1];
sx q[1];
rz(-2.5187056) q[1];
x q[2];
rz(1.5434389) q[3];
sx q[3];
rz(-1.1809071) q[3];
sx q[3];
rz(1.1297806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.61218843) q[2];
sx q[2];
rz(-0.89451423) q[2];
sx q[2];
rz(0.55076304) q[2];
rz(-2.4380056) q[3];
sx q[3];
rz(-1.0102605) q[3];
sx q[3];
rz(-1.5065058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7228912) q[0];
sx q[0];
rz(-2.7846865) q[0];
sx q[0];
rz(3.0974467) q[0];
rz(1.528953) q[1];
sx q[1];
rz(-1.2395369) q[1];
sx q[1];
rz(-2.3619161) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3846489) q[0];
sx q[0];
rz(-2.1063519) q[0];
sx q[0];
rz(1.4063565) q[0];
rz(-pi) q[1];
x q[1];
rz(0.81515628) q[2];
sx q[2];
rz(-1.2218468) q[2];
sx q[2];
rz(-2.7100035) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.30291468) q[1];
sx q[1];
rz(-1.0524787) q[1];
sx q[1];
rz(1.7811437) q[1];
rz(0.78124222) q[3];
sx q[3];
rz(-1.9133854) q[3];
sx q[3];
rz(0.70171802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.49446517) q[2];
sx q[2];
rz(-0.72138849) q[2];
sx q[2];
rz(-1.5987827) q[2];
rz(-0.70458448) q[3];
sx q[3];
rz(-1.5141809) q[3];
sx q[3];
rz(-0.012044756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.71173944) q[0];
sx q[0];
rz(-1.5562417) q[0];
sx q[0];
rz(-0.65162311) q[0];
rz(-1.7977057) q[1];
sx q[1];
rz(-1.6603036) q[1];
sx q[1];
rz(2.4659326) q[1];
rz(-2.007953) q[2];
sx q[2];
rz(-1.2883678) q[2];
sx q[2];
rz(-1.012445) q[2];
rz(0.6380973) q[3];
sx q[3];
rz(-1.4641855) q[3];
sx q[3];
rz(-2.2104213) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
