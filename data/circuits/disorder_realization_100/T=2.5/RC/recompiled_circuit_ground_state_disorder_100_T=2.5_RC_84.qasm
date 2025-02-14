OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7772726) q[0];
sx q[0];
rz(-0.84982818) q[0];
sx q[0];
rz(1.6027236) q[0];
rz(-0.76378769) q[1];
sx q[1];
rz(-0.053215947) q[1];
sx q[1];
rz(-2.3071212) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.065782308) q[0];
sx q[0];
rz(-0.47496034) q[0];
sx q[0];
rz(1.4592828) q[0];
rz(0.75196377) q[2];
sx q[2];
rz(-2.3533216) q[2];
sx q[2];
rz(-1.0984818) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.149772) q[1];
sx q[1];
rz(-1.705994) q[1];
sx q[1];
rz(-0.97414952) q[1];
rz(0.90596195) q[3];
sx q[3];
rz(-2.6506781) q[3];
sx q[3];
rz(-0.57074805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6795464) q[2];
sx q[2];
rz(-1.5908238) q[2];
sx q[2];
rz(-0.66590434) q[2];
rz(-0.42689651) q[3];
sx q[3];
rz(-0.96564966) q[3];
sx q[3];
rz(-2.9396333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17215984) q[0];
sx q[0];
rz(-0.96123022) q[0];
sx q[0];
rz(2.6184671) q[0];
rz(-2.7223749) q[1];
sx q[1];
rz(-0.20877498) q[1];
sx q[1];
rz(2.8255393) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5764874) q[0];
sx q[0];
rz(-2.2448328) q[0];
sx q[0];
rz(1.7934858) q[0];
x q[1];
rz(-0.14618539) q[2];
sx q[2];
rz(-1.5954895) q[2];
sx q[2];
rz(3.0865412) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.167724) q[1];
sx q[1];
rz(-1.8710724) q[1];
sx q[1];
rz(2.1435277) q[1];
rz(-0.42556588) q[3];
sx q[3];
rz(-1.4168973) q[3];
sx q[3];
rz(-0.16395149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.54059625) q[2];
sx q[2];
rz(-2.5935846) q[2];
sx q[2];
rz(-2.0151286) q[2];
rz(2.2362041) q[3];
sx q[3];
rz(-1.0631881) q[3];
sx q[3];
rz(1.6006914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0970522) q[0];
sx q[0];
rz(-1.4690761) q[0];
sx q[0];
rz(1.4918406) q[0];
rz(-0.39332321) q[1];
sx q[1];
rz(-1.2231188) q[1];
sx q[1];
rz(-2.9755859) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9764023) q[0];
sx q[0];
rz(-0.87346625) q[0];
sx q[0];
rz(1.9843452) q[0];
x q[1];
rz(-2.2691378) q[2];
sx q[2];
rz(-1.7042245) q[2];
sx q[2];
rz(-2.7492571) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5854726) q[1];
sx q[1];
rz(-0.86143805) q[1];
sx q[1];
rz(1.8147537) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7083206) q[3];
sx q[3];
rz(-1.6973933) q[3];
sx q[3];
rz(1.8357852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4094706) q[2];
sx q[2];
rz(-0.25264838) q[2];
sx q[2];
rz(-2.9518993) q[2];
rz(-0.41483748) q[3];
sx q[3];
rz(-1.8388351) q[3];
sx q[3];
rz(-3.0366669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46284997) q[0];
sx q[0];
rz(-2.3609556) q[0];
sx q[0];
rz(-1.2586052) q[0];
rz(1.3858093) q[1];
sx q[1];
rz(-0.36711991) q[1];
sx q[1];
rz(-0.78305903) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.258062) q[0];
sx q[0];
rz(-1.5386762) q[0];
sx q[0];
rz(1.5954992) q[0];
rz(-pi) q[1];
rz(0.80863805) q[2];
sx q[2];
rz(-1.5735192) q[2];
sx q[2];
rz(2.115094) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.051797051) q[1];
sx q[1];
rz(-2.0690114) q[1];
sx q[1];
rz(2.5913622) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5989877) q[3];
sx q[3];
rz(-1.2148093) q[3];
sx q[3];
rz(0.062075768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.9328402) q[2];
sx q[2];
rz(-1.1609266) q[2];
sx q[2];
rz(1.0932385) q[2];
rz(-0.42539635) q[3];
sx q[3];
rz(-1.2109816) q[3];
sx q[3];
rz(1.0824599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071674034) q[0];
sx q[0];
rz(-2.2921083) q[0];
sx q[0];
rz(1.7639532) q[0];
rz(0.90179521) q[1];
sx q[1];
rz(-1.8865562) q[1];
sx q[1];
rz(-0.76595438) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1655584) q[0];
sx q[0];
rz(-1.9355825) q[0];
sx q[0];
rz(-2.4976612) q[0];
rz(-pi) q[1];
rz(-1.2754006) q[2];
sx q[2];
rz(-0.89933853) q[2];
sx q[2];
rz(1.7341933) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.9498038) q[1];
sx q[1];
rz(-1.7226137) q[1];
sx q[1];
rz(-2.6372974) q[1];
rz(-pi) q[2];
rz(-2.125191) q[3];
sx q[3];
rz(-0.50339744) q[3];
sx q[3];
rz(2.7254846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.19096601) q[2];
sx q[2];
rz(-2.2910255) q[2];
sx q[2];
rz(-1.4858049) q[2];
rz(2.0048257) q[3];
sx q[3];
rz(-2.6906689) q[3];
sx q[3];
rz(2.3740812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.23761) q[0];
sx q[0];
rz(-0.5466277) q[0];
sx q[0];
rz(-2.9055063) q[0];
rz(1.4769332) q[1];
sx q[1];
rz(-1.4918985) q[1];
sx q[1];
rz(-1.1621071) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2598545) q[0];
sx q[0];
rz(-0.82603541) q[0];
sx q[0];
rz(-1.8709454) q[0];
rz(-pi) q[1];
rz(-2.7323551) q[2];
sx q[2];
rz(-1.6820568) q[2];
sx q[2];
rz(-1.2986599) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9182738) q[1];
sx q[1];
rz(-1.2122723) q[1];
sx q[1];
rz(0.29488342) q[1];
rz(0.052836494) q[3];
sx q[3];
rz(-1.5457258) q[3];
sx q[3];
rz(0.79475886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.37940019) q[2];
sx q[2];
rz(-0.47319943) q[2];
sx q[2];
rz(0.0089448373) q[2];
rz(1.829932) q[3];
sx q[3];
rz(-2.0904171) q[3];
sx q[3];
rz(-1.8478307) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0084429) q[0];
sx q[0];
rz(-2.382998) q[0];
sx q[0];
rz(-1.5593418) q[0];
rz(-0.82296952) q[1];
sx q[1];
rz(-2.6428849) q[1];
sx q[1];
rz(1.4272326) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2555851) q[0];
sx q[0];
rz(-0.56047165) q[0];
sx q[0];
rz(0.94985243) q[0];
x q[1];
rz(0.64325579) q[2];
sx q[2];
rz(-1.3998886) q[2];
sx q[2];
rz(-2.3115445) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.848915) q[1];
sx q[1];
rz(-2.1452322) q[1];
sx q[1];
rz(0.019197063) q[1];
rz(-0.41540878) q[3];
sx q[3];
rz(-1.4671638) q[3];
sx q[3];
rz(-1.2670011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3048749) q[2];
sx q[2];
rz(-2.0676401) q[2];
sx q[2];
rz(2.7433336) q[2];
rz(2.7598329) q[3];
sx q[3];
rz(-1.441381) q[3];
sx q[3];
rz(2.1329342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1036103) q[0];
sx q[0];
rz(-1.7033662) q[0];
sx q[0];
rz(1.9110699) q[0];
rz(1.1331753) q[1];
sx q[1];
rz(-0.76825348) q[1];
sx q[1];
rz(0.59655985) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2418848) q[0];
sx q[0];
rz(-1.5615789) q[0];
sx q[0];
rz(-2.8901744) q[0];
rz(-1.9279154) q[2];
sx q[2];
rz(-2.6304863) q[2];
sx q[2];
rz(-0.36506436) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.34932) q[1];
sx q[1];
rz(-1.8304003) q[1];
sx q[1];
rz(-0.97430299) q[1];
x q[2];
rz(1.1876218) q[3];
sx q[3];
rz(-2.1541693) q[3];
sx q[3];
rz(0.72900984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.29264984) q[2];
sx q[2];
rz(-0.611648) q[2];
sx q[2];
rz(-1.4673648) q[2];
rz(-0.26662207) q[3];
sx q[3];
rz(-0.98054236) q[3];
sx q[3];
rz(0.22661859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.4416696) q[0];
sx q[0];
rz(-3.141098) q[0];
sx q[0];
rz(1.1215522) q[0];
rz(-0.17555217) q[1];
sx q[1];
rz(-0.67070621) q[1];
sx q[1];
rz(-0.49351722) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7132015) q[0];
sx q[0];
rz(-1.6294217) q[0];
sx q[0];
rz(-1.4992264) q[0];
x q[1];
rz(1.1404806) q[2];
sx q[2];
rz(-2.3386933) q[2];
sx q[2];
rz(-1.0969298) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.15010506) q[1];
sx q[1];
rz(-0.59044033) q[1];
sx q[1];
rz(-1.9976678) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0974985) q[3];
sx q[3];
rz(-1.4933503) q[3];
sx q[3];
rz(0.29260422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.4404122) q[2];
sx q[2];
rz(-2.7737325) q[2];
sx q[2];
rz(-0.58604798) q[2];
rz(-1.690257) q[3];
sx q[3];
rz(-1.3380932) q[3];
sx q[3];
rz(2.6293248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.0217836) q[0];
sx q[0];
rz(-1.021294) q[0];
sx q[0];
rz(2.9721066) q[0];
rz(-2.4985864) q[1];
sx q[1];
rz(-2.2702859) q[1];
sx q[1];
rz(-2.6298584) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.046233072) q[0];
sx q[0];
rz(-1.4698782) q[0];
sx q[0];
rz(-1.4907229) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8795661) q[2];
sx q[2];
rz(-1.6341156) q[2];
sx q[2];
rz(2.2514907) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4525254) q[1];
sx q[1];
rz(-0.91740174) q[1];
sx q[1];
rz(2.9513068) q[1];
x q[2];
rz(-1.9668642) q[3];
sx q[3];
rz(-0.549256) q[3];
sx q[3];
rz(1.6793457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.52342755) q[2];
sx q[2];
rz(-0.7343505) q[2];
sx q[2];
rz(0.0283453) q[2];
rz(0.47447765) q[3];
sx q[3];
rz(-0.35895434) q[3];
sx q[3];
rz(2.1019782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9356553) q[0];
sx q[0];
rz(-1.7556998) q[0];
sx q[0];
rz(2.6082267) q[0];
rz(-0.43279303) q[1];
sx q[1];
rz(-1.5569893) q[1];
sx q[1];
rz(2.5098262) q[1];
rz(-2.8001187) q[2];
sx q[2];
rz(-0.22814422) q[2];
sx q[2];
rz(3.0244915) q[2];
rz(-1.9267052) q[3];
sx q[3];
rz(-0.22161814) q[3];
sx q[3];
rz(1.591481) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
