OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7423695) q[0];
sx q[0];
rz(2.7288781) q[0];
sx q[0];
rz(5.4469845) q[0];
rz(-2.3669481) q[1];
sx q[1];
rz(6.3451938) q[1];
sx q[1];
rz(11.558029) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8275534) q[0];
sx q[0];
rz(-1.9441105) q[0];
sx q[0];
rz(-0.084777431) q[0];
x q[1];
rz(0.37990976) q[2];
sx q[2];
rz(-1.3515944) q[2];
sx q[2];
rz(2.4394112) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5656149) q[1];
sx q[1];
rz(-1.4139785) q[1];
sx q[1];
rz(-0.48139907) q[1];
rz(0.020691607) q[3];
sx q[3];
rz(-1.3776099) q[3];
sx q[3];
rz(2.8789946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8334373) q[2];
sx q[2];
rz(-0.94352949) q[2];
sx q[2];
rz(-0.65307871) q[2];
rz(0.11450442) q[3];
sx q[3];
rz(-1.7479892) q[3];
sx q[3];
rz(1.0623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7750875) q[0];
sx q[0];
rz(-1.0307642) q[0];
sx q[0];
rz(1.7979701) q[0];
rz(1.1603181) q[1];
sx q[1];
rz(-0.41028816) q[1];
sx q[1];
rz(0.62058273) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68408332) q[0];
sx q[0];
rz(-2.7433222) q[0];
sx q[0];
rz(-1.0578367) q[0];
rz(-pi) q[1];
rz(1.1158022) q[2];
sx q[2];
rz(-2.4206851) q[2];
sx q[2];
rz(1.4269331) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0936135) q[1];
sx q[1];
rz(-1.9318214) q[1];
sx q[1];
rz(0.71012902) q[1];
rz(-pi) q[2];
rz(-2.85895) q[3];
sx q[3];
rz(-0.97393721) q[3];
sx q[3];
rz(-2.7925036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4429984) q[2];
sx q[2];
rz(-1.6220762) q[2];
sx q[2];
rz(-0.2300187) q[2];
rz(1.5551785) q[3];
sx q[3];
rz(-1.9929726) q[3];
sx q[3];
rz(-2.8659099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8860633) q[0];
sx q[0];
rz(-2.7524188) q[0];
sx q[0];
rz(1.080876) q[0];
rz(-2.4983662) q[1];
sx q[1];
rz(-0.79935646) q[1];
sx q[1];
rz(3.0562775) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9935308) q[0];
sx q[0];
rz(-2.2087653) q[0];
sx q[0];
rz(-0.83420269) q[0];
rz(0.32344748) q[2];
sx q[2];
rz(-2.5451042) q[2];
sx q[2];
rz(1.310854) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4115632) q[1];
sx q[1];
rz(-1.7021028) q[1];
sx q[1];
rz(0.16741411) q[1];
rz(-pi) q[2];
rz(-1.6162698) q[3];
sx q[3];
rz(-0.95004987) q[3];
sx q[3];
rz(2.5705702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9366511) q[2];
sx q[2];
rz(-3.1320269) q[2];
sx q[2];
rz(-0.19634518) q[2];
rz(-0.28228545) q[3];
sx q[3];
rz(-2.0245656) q[3];
sx q[3];
rz(-2.096368) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1333756) q[0];
sx q[0];
rz(-0.5468002) q[0];
sx q[0];
rz(-2.7561482) q[0];
rz(1.7281945) q[1];
sx q[1];
rz(-1.9368659) q[1];
sx q[1];
rz(-2.8916496) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2521533) q[0];
sx q[0];
rz(-1.982054) q[0];
sx q[0];
rz(-0.90775872) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6310299) q[2];
sx q[2];
rz(-1.4830198) q[2];
sx q[2];
rz(-2.6812832) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.39482597) q[1];
sx q[1];
rz(-1.4020846) q[1];
sx q[1];
rz(-0.49612237) q[1];
rz(-pi) q[2];
rz(-0.79367359) q[3];
sx q[3];
rz(-1.4303615) q[3];
sx q[3];
rz(0.92257231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.46127737) q[2];
sx q[2];
rz(-1.8467434) q[2];
sx q[2];
rz(0.2363905) q[2];
rz(1.8810898) q[3];
sx q[3];
rz(-2.8915296) q[3];
sx q[3];
rz(-0.67729706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8054955) q[0];
sx q[0];
rz(-1.6039055) q[0];
sx q[0];
rz(0.48511037) q[0];
rz(1.9612034) q[1];
sx q[1];
rz(-1.5472658) q[1];
sx q[1];
rz(-0.83650437) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77959767) q[0];
sx q[0];
rz(-1.2149724) q[0];
sx q[0];
rz(-1.9076288) q[0];
x q[1];
rz(-2.8233158) q[2];
sx q[2];
rz(-1.1851382) q[2];
sx q[2];
rz(1.1952656) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6753458) q[1];
sx q[1];
rz(-2.2398178) q[1];
sx q[1];
rz(-0.061208486) q[1];
rz(-pi) q[2];
rz(-0.4007217) q[3];
sx q[3];
rz(-1.5127137) q[3];
sx q[3];
rz(0.42089128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.89852077) q[2];
sx q[2];
rz(-2.409755) q[2];
sx q[2];
rz(-2.3720429) q[2];
rz(0.92929333) q[3];
sx q[3];
rz(-1.1309036) q[3];
sx q[3];
rz(0.23269674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-1.1997851) q[0];
sx q[0];
rz(-0.48518825) q[0];
sx q[0];
rz(0.76882452) q[0];
rz(-1.3517514) q[1];
sx q[1];
rz(-0.47305802) q[1];
sx q[1];
rz(2.2519055) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6883967) q[0];
sx q[0];
rz(-1.9510117) q[0];
sx q[0];
rz(0.225747) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5403929) q[2];
sx q[2];
rz(-0.044375751) q[2];
sx q[2];
rz(2.6868683) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5701712) q[1];
sx q[1];
rz(-1.716905) q[1];
sx q[1];
rz(-3.094638) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4558639) q[3];
sx q[3];
rz(-1.1224261) q[3];
sx q[3];
rz(-1.861972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6512904) q[2];
sx q[2];
rz(-1.7040665) q[2];
sx q[2];
rz(1.5865405) q[2];
rz(1.2706612) q[3];
sx q[3];
rz(-2.7226166) q[3];
sx q[3];
rz(0.13203013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9689869) q[0];
sx q[0];
rz(-2.0033328) q[0];
sx q[0];
rz(-1.1440811) q[0];
rz(0.83089685) q[1];
sx q[1];
rz(-1.7832719) q[1];
sx q[1];
rz(-2.3544618) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4517196) q[0];
sx q[0];
rz(-1.7053118) q[0];
sx q[0];
rz(-3.0089507) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.5412124) q[2];
sx q[2];
rz(-2.1809309) q[2];
sx q[2];
rz(1.0224563) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5435782) q[1];
sx q[1];
rz(-2.0435963) q[1];
sx q[1];
rz(0.531457) q[1];
rz(-0.55946405) q[3];
sx q[3];
rz(-1.2351523) q[3];
sx q[3];
rz(-1.0192724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1181011) q[2];
sx q[2];
rz(-2.1716437) q[2];
sx q[2];
rz(-0.033163158) q[2];
rz(-1.5983332) q[3];
sx q[3];
rz(-0.71591806) q[3];
sx q[3];
rz(-1.5020812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.892266) q[0];
sx q[0];
rz(-1.5557657) q[0];
sx q[0];
rz(1.7484885) q[0];
rz(2.6847367) q[1];
sx q[1];
rz(-1.9970048) q[1];
sx q[1];
rz(2.0410062) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4286713) q[0];
sx q[0];
rz(-1.5690737) q[0];
sx q[0];
rz(0.07600204) q[0];
rz(0.92488168) q[2];
sx q[2];
rz(-0.98306984) q[2];
sx q[2];
rz(-2.4921668) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9683084) q[1];
sx q[1];
rz(-1.5368764) q[1];
sx q[1];
rz(0.90356253) q[1];
rz(-pi) q[2];
rz(-1.4310775) q[3];
sx q[3];
rz(-1.1828239) q[3];
sx q[3];
rz(-3.1041077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0988203) q[2];
sx q[2];
rz(-2.6440812) q[2];
sx q[2];
rz(-1.5808606) q[2];
rz(2.3780195) q[3];
sx q[3];
rz(-1.5127425) q[3];
sx q[3];
rz(-1.0360576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40912691) q[0];
sx q[0];
rz(-2.3865073) q[0];
sx q[0];
rz(-2.8501046) q[0];
rz(1.905722) q[1];
sx q[1];
rz(-1.9449077) q[1];
sx q[1];
rz(-0.12289563) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6345469) q[0];
sx q[0];
rz(-1.5193918) q[0];
sx q[0];
rz(2.794751) q[0];
x q[1];
rz(-2.2332623) q[2];
sx q[2];
rz(-0.35229063) q[2];
sx q[2];
rz(0.34991821) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.025561573) q[1];
sx q[1];
rz(-1.681038) q[1];
sx q[1];
rz(-1.362839) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4261774) q[3];
sx q[3];
rz(-1.5767617) q[3];
sx q[3];
rz(2.1534513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5128532) q[2];
sx q[2];
rz(-1.3924007) q[2];
sx q[2];
rz(1.0635618) q[2];
rz(-2.9641446) q[3];
sx q[3];
rz(-0.95824233) q[3];
sx q[3];
rz(-2.1232846) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7040831) q[0];
sx q[0];
rz(-2.6884485) q[0];
sx q[0];
rz(-1.4310687) q[0];
rz(-0.60072947) q[1];
sx q[1];
rz(-0.44980106) q[1];
sx q[1];
rz(-0.61753714) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3356614) q[0];
sx q[0];
rz(-1.7863723) q[0];
sx q[0];
rz(3.0144601) q[0];
rz(-pi) q[1];
rz(0.82771684) q[2];
sx q[2];
rz(-1.9736145) q[2];
sx q[2];
rz(-3.0981491) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.75467089) q[1];
sx q[1];
rz(-0.53334177) q[1];
sx q[1];
rz(-1.6120595) q[1];
rz(2.6838949) q[3];
sx q[3];
rz(-2.5182596) q[3];
sx q[3];
rz(1.8231572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.90126976) q[2];
sx q[2];
rz(-2.1412886) q[2];
sx q[2];
rz(-2.1642302) q[2];
rz(2.0591002) q[3];
sx q[3];
rz(-2.2668224) q[3];
sx q[3];
rz(3.0860743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8025773) q[0];
sx q[0];
rz(-0.76887283) q[0];
sx q[0];
rz(-0.73706891) q[0];
rz(3.0958685) q[1];
sx q[1];
rz(-2.3241691) q[1];
sx q[1];
rz(-1.587422) q[1];
rz(-2.634766) q[2];
sx q[2];
rz(-2.0438556) q[2];
sx q[2];
rz(0.023509916) q[2];
rz(0.65683881) q[3];
sx q[3];
rz(-1.7056864) q[3];
sx q[3];
rz(1.0581072) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
