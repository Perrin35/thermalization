OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4274974) q[0];
sx q[0];
rz(-0.56199718) q[0];
sx q[0];
rz(-2.9105817) q[0];
rz(-2.8958939) q[1];
sx q[1];
rz(-2.6872771) q[1];
sx q[1];
rz(1.8543724) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2184233) q[0];
sx q[0];
rz(-1.2997264) q[0];
sx q[0];
rz(-2.0725033) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76164772) q[2];
sx q[2];
rz(-0.33312329) q[2];
sx q[2];
rz(-1.3093349) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2778863) q[1];
sx q[1];
rz(-1.4046298) q[1];
sx q[1];
rz(-1.4708797) q[1];
rz(1.1679959) q[3];
sx q[3];
rz(-2.9514545) q[3];
sx q[3];
rz(-2.8791752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7370558) q[2];
sx q[2];
rz(-0.70713592) q[2];
sx q[2];
rz(2.7810968) q[2];
rz(-1.1245842) q[3];
sx q[3];
rz(-2.093061) q[3];
sx q[3];
rz(1.2688961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8812113) q[0];
sx q[0];
rz(-0.17558782) q[0];
sx q[0];
rz(0.65688175) q[0];
rz(0.33292133) q[1];
sx q[1];
rz(-1.0924783) q[1];
sx q[1];
rz(-0.93516707) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012197709) q[0];
sx q[0];
rz(-3.0327065) q[0];
sx q[0];
rz(-0.28277855) q[0];
rz(-pi) q[1];
x q[1];
rz(0.25721154) q[2];
sx q[2];
rz(-1.6597444) q[2];
sx q[2];
rz(-2.6572029) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.540024) q[1];
sx q[1];
rz(-1.4238796) q[1];
sx q[1];
rz(1.6107035) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3627429) q[3];
sx q[3];
rz(-1.2223772) q[3];
sx q[3];
rz(1.5294242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8344581) q[2];
sx q[2];
rz(-1.6259401) q[2];
sx q[2];
rz(1.2163986) q[2];
rz(0.52792102) q[3];
sx q[3];
rz(-1.0390176) q[3];
sx q[3];
rz(1.2549887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6127748) q[0];
sx q[0];
rz(-2.3154494) q[0];
sx q[0];
rz(-0.58498996) q[0];
rz(3.0168369) q[1];
sx q[1];
rz(-2.545732) q[1];
sx q[1];
rz(-1.4488719) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0199658) q[0];
sx q[0];
rz(-0.0040071132) q[0];
sx q[0];
rz(0.83929707) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.072183822) q[2];
sx q[2];
rz(-1.5572786) q[2];
sx q[2];
rz(1.9428321) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3347149) q[1];
sx q[1];
rz(-0.92744614) q[1];
sx q[1];
rz(-0.54912864) q[1];
rz(-pi) q[2];
rz(2.8040941) q[3];
sx q[3];
rz(-1.322116) q[3];
sx q[3];
rz(1.6832222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2584194) q[2];
sx q[2];
rz(-2.1397739) q[2];
sx q[2];
rz(-1.4460571) q[2];
rz(0.7575194) q[3];
sx q[3];
rz(-1.5599374) q[3];
sx q[3];
rz(-0.83938804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.7739173) q[0];
sx q[0];
rz(-0.92905074) q[0];
sx q[0];
rz(1.0876592) q[0];
rz(-2.7663973) q[1];
sx q[1];
rz(-2.0162069) q[1];
sx q[1];
rz(2.1813724) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4488569) q[0];
sx q[0];
rz(-1.6105284) q[0];
sx q[0];
rz(1.8360774) q[0];
rz(-pi) q[1];
rz(1.8619991) q[2];
sx q[2];
rz(-3.0377977) q[2];
sx q[2];
rz(-1.7218931) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6672047) q[1];
sx q[1];
rz(-1.9992113) q[1];
sx q[1];
rz(-1.0775671) q[1];
x q[2];
rz(2.8720565) q[3];
sx q[3];
rz(-1.4785267) q[3];
sx q[3];
rz(2.5470981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.20597657) q[2];
sx q[2];
rz(-1.2674067) q[2];
sx q[2];
rz(1.6440294) q[2];
rz(2.2802672) q[3];
sx q[3];
rz(-2.438811) q[3];
sx q[3];
rz(0.45708814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.977026) q[0];
sx q[0];
rz(-2.470546) q[0];
sx q[0];
rz(0.60428756) q[0];
rz(-1.1445649) q[1];
sx q[1];
rz(-1.4860169) q[1];
sx q[1];
rz(-2.7191275) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2328932) q[0];
sx q[0];
rz(-0.18747231) q[0];
sx q[0];
rz(1.1195681) q[0];
rz(-0.79549148) q[2];
sx q[2];
rz(-2.4078712) q[2];
sx q[2];
rz(-0.41662859) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8026655) q[1];
sx q[1];
rz(-2.4484854) q[1];
sx q[1];
rz(3.0285804) q[1];
x q[2];
rz(1.2286387) q[3];
sx q[3];
rz(-1.9316202) q[3];
sx q[3];
rz(-2.7865041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3182688) q[2];
sx q[2];
rz(-2.4906929) q[2];
sx q[2];
rz(0.65625119) q[2];
rz(1.5308135) q[3];
sx q[3];
rz(-1.2698413) q[3];
sx q[3];
rz(-2.5608565) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4055279) q[0];
sx q[0];
rz(-1.0117714) q[0];
sx q[0];
rz(2.5166125) q[0];
rz(0.75195733) q[1];
sx q[1];
rz(-1.0686921) q[1];
sx q[1];
rz(1.5350852) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2921755) q[0];
sx q[0];
rz(-1.8596974) q[0];
sx q[0];
rz(-1.6221415) q[0];
x q[1];
rz(-2.4703896) q[2];
sx q[2];
rz(-0.83121383) q[2];
sx q[2];
rz(2.9321456) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5175486) q[1];
sx q[1];
rz(-2.4772948) q[1];
sx q[1];
rz(-2.7655927) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5967136) q[3];
sx q[3];
rz(-2.0996465) q[3];
sx q[3];
rz(1.8008055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.1890761) q[2];
sx q[2];
rz(-1.357888) q[2];
sx q[2];
rz(0.9291741) q[2];
rz(-2.3164228) q[3];
sx q[3];
rz(-1.5114096) q[3];
sx q[3];
rz(1.1024124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5850942) q[0];
sx q[0];
rz(-2.3346021) q[0];
sx q[0];
rz(-1.7671385) q[0];
rz(0.51042405) q[1];
sx q[1];
rz(-2.3511032) q[1];
sx q[1];
rz(0.62320954) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40996288) q[0];
sx q[0];
rz(-1.7222758) q[0];
sx q[0];
rz(1.234647) q[0];
rz(0.34131949) q[2];
sx q[2];
rz(-0.53487294) q[2];
sx q[2];
rz(-2.24077) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0125) q[1];
sx q[1];
rz(-2.9575037) q[1];
sx q[1];
rz(-0.91839183) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74759746) q[3];
sx q[3];
rz(-1.4211402) q[3];
sx q[3];
rz(-0.29779321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0906543) q[2];
sx q[2];
rz(-2.3263558) q[2];
sx q[2];
rz(1.9742924) q[2];
rz(0.022631571) q[3];
sx q[3];
rz(-0.52112094) q[3];
sx q[3];
rz(-1.1289271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.8959494) q[0];
sx q[0];
rz(-1.1351981) q[0];
sx q[0];
rz(2.1242712) q[0];
rz(1.3941049) q[1];
sx q[1];
rz(-2.930495) q[1];
sx q[1];
rz(0.62754935) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2557954) q[0];
sx q[0];
rz(-1.9667407) q[0];
sx q[0];
rz(-0.93314472) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4862138) q[2];
sx q[2];
rz(-2.248317) q[2];
sx q[2];
rz(-1.5280452) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7146475) q[1];
sx q[1];
rz(-1.8083739) q[1];
sx q[1];
rz(-1.7725043) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5379337) q[3];
sx q[3];
rz(-2.5839879) q[3];
sx q[3];
rz(1.738036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.62577406) q[2];
sx q[2];
rz(-2.5067582) q[2];
sx q[2];
rz(0.41075692) q[2];
rz(3.0913894) q[3];
sx q[3];
rz(-1.8886731) q[3];
sx q[3];
rz(-3.0661809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5589767) q[0];
sx q[0];
rz(-0.89685431) q[0];
sx q[0];
rz(-2.8726752) q[0];
rz(0.31002632) q[1];
sx q[1];
rz(-1.0870442) q[1];
sx q[1];
rz(0.1300098) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.01154218) q[0];
sx q[0];
rz(-0.39549144) q[0];
sx q[0];
rz(-2.1756835) q[0];
rz(-2.347888) q[2];
sx q[2];
rz(-2.6598601) q[2];
sx q[2];
rz(-0.46221737) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8925608) q[1];
sx q[1];
rz(-1.4098123) q[1];
sx q[1];
rz(-0.56610961) q[1];
rz(0.49874108) q[3];
sx q[3];
rz(-2.1864656) q[3];
sx q[3];
rz(-1.2241253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2076063) q[2];
sx q[2];
rz(-2.3748368) q[2];
sx q[2];
rz(3.0450191) q[2];
rz(1.5038331) q[3];
sx q[3];
rz(-1.8042754) q[3];
sx q[3];
rz(-0.79969978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0132975) q[0];
sx q[0];
rz(-0.18823637) q[0];
sx q[0];
rz(-0.53303322) q[0];
rz(-0.036272613) q[1];
sx q[1];
rz(-2.3590922) q[1];
sx q[1];
rz(-1.689555) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5851327) q[0];
sx q[0];
rz(-1.0074573) q[0];
sx q[0];
rz(-2.6908532) q[0];
rz(-0.31970892) q[2];
sx q[2];
rz(-0.33505782) q[2];
sx q[2];
rz(2.8220995) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2382792) q[1];
sx q[1];
rz(-2.7809445) q[1];
sx q[1];
rz(-3.0566932) q[1];
rz(-pi) q[2];
x q[2];
rz(2.96569) q[3];
sx q[3];
rz(-0.47187343) q[3];
sx q[3];
rz(-0.37341082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4286246) q[2];
sx q[2];
rz(-0.73178256) q[2];
sx q[2];
rz(3.1089605) q[2];
rz(0.84515682) q[3];
sx q[3];
rz(-1.2266351) q[3];
sx q[3];
rz(-2.0613861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7964771) q[0];
sx q[0];
rz(-0.65237541) q[0];
sx q[0];
rz(2.8271578) q[0];
rz(2.3445917) q[1];
sx q[1];
rz(-2.2140257) q[1];
sx q[1];
rz(0.066233403) q[1];
rz(-2.6673139) q[2];
sx q[2];
rz(-0.58816345) q[2];
sx q[2];
rz(0.11810398) q[2];
rz(-1.4440404) q[3];
sx q[3];
rz(-0.76382617) q[3];
sx q[3];
rz(2.7504117) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
