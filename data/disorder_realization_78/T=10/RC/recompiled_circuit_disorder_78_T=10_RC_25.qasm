OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.55943638) q[0];
sx q[0];
rz(3.732582) q[0];
sx q[0];
rz(8.8413722) q[0];
rz(-0.18435873) q[1];
sx q[1];
rz(-2.1579722) q[1];
sx q[1];
rz(0.89259994) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28461449) q[0];
sx q[0];
rz(-2.1224788) q[0];
sx q[0];
rz(2.7906228) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66944389) q[2];
sx q[2];
rz(-1.9813683) q[2];
sx q[2];
rz(0.66561156) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0639227) q[1];
sx q[1];
rz(-2.3706672) q[1];
sx q[1];
rz(-2.0248807) q[1];
rz(-pi) q[2];
rz(1.5276018) q[3];
sx q[3];
rz(-2.3213263) q[3];
sx q[3];
rz(0.84212069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2261752) q[2];
sx q[2];
rz(-1.5390652) q[2];
sx q[2];
rz(1.9809451) q[2];
rz(-2.9246269) q[3];
sx q[3];
rz(-2.6187077) q[3];
sx q[3];
rz(-1.0552361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.083667) q[0];
sx q[0];
rz(-2.1766429) q[0];
sx q[0];
rz(-0.57587409) q[0];
rz(1.8946164) q[1];
sx q[1];
rz(-1.2966825) q[1];
sx q[1];
rz(1.974568) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8407699) q[0];
sx q[0];
rz(-1.5123899) q[0];
sx q[0];
rz(1.3129243) q[0];
rz(-pi) q[1];
rz(-2.1115233) q[2];
sx q[2];
rz(-1.3582555) q[2];
sx q[2];
rz(-2.9589047) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.45707073) q[1];
sx q[1];
rz(-1.6324537) q[1];
sx q[1];
rz(-2.3350299) q[1];
x q[2];
rz(-2.0608276) q[3];
sx q[3];
rz(-0.96066517) q[3];
sx q[3];
rz(-0.78435635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1077659) q[2];
sx q[2];
rz(-1.9627389) q[2];
sx q[2];
rz(-2.9272184) q[2];
rz(3.068148) q[3];
sx q[3];
rz(-0.44973222) q[3];
sx q[3];
rz(-0.28081056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.4784933) q[0];
sx q[0];
rz(-0.61613023) q[0];
sx q[0];
rz(1.2269155) q[0];
rz(-0.40027174) q[1];
sx q[1];
rz(-1.8881256) q[1];
sx q[1];
rz(-2.1267557) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6944511) q[0];
sx q[0];
rz(-1.5146291) q[0];
sx q[0];
rz(2.4980314) q[0];
x q[1];
rz(-2.8995908) q[2];
sx q[2];
rz(-1.9654462) q[2];
sx q[2];
rz(1.637527) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9008357) q[1];
sx q[1];
rz(-1.6512617) q[1];
sx q[1];
rz(-2.1686173) q[1];
x q[2];
rz(-0.26168163) q[3];
sx q[3];
rz(-0.90173429) q[3];
sx q[3];
rz(-0.79479587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0959452) q[2];
sx q[2];
rz(-2.084705) q[2];
sx q[2];
rz(0.8992368) q[2];
rz(2.4441161) q[3];
sx q[3];
rz(-1.8557502) q[3];
sx q[3];
rz(-2.5337059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4863481) q[0];
sx q[0];
rz(-1.0826033) q[0];
sx q[0];
rz(-2.7096601) q[0];
rz(0.63255429) q[1];
sx q[1];
rz(-2.7245941) q[1];
sx q[1];
rz(-0.63582173) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0974554) q[0];
sx q[0];
rz(-2.6090528) q[0];
sx q[0];
rz(-1.9299279) q[0];
rz(-pi) q[1];
rz(2.9767838) q[2];
sx q[2];
rz(-0.85068446) q[2];
sx q[2];
rz(2.1317496) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1078474) q[1];
sx q[1];
rz(-0.81106942) q[1];
sx q[1];
rz(-1.8268405) q[1];
rz(2.8549529) q[3];
sx q[3];
rz(-1.7754103) q[3];
sx q[3];
rz(-0.74433792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2137961) q[2];
sx q[2];
rz(-2.9979604) q[2];
sx q[2];
rz(-2.4528743) q[2];
rz(0.33411807) q[3];
sx q[3];
rz(-1.170661) q[3];
sx q[3];
rz(-2.9978602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14389811) q[0];
sx q[0];
rz(-0.69902885) q[0];
sx q[0];
rz(-2.0671663) q[0];
rz(-2.396446) q[1];
sx q[1];
rz(-1.6586168) q[1];
sx q[1];
rz(0.27854663) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4771839) q[0];
sx q[0];
rz(-1.1468977) q[0];
sx q[0];
rz(0.82188481) q[0];
rz(-pi) q[1];
rz(-1.6804382) q[2];
sx q[2];
rz(-0.96193681) q[2];
sx q[2];
rz(1.7313752) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6812233) q[1];
sx q[1];
rz(-1.9378098) q[1];
sx q[1];
rz(0.91194921) q[1];
rz(-pi) q[2];
rz(0.99028011) q[3];
sx q[3];
rz(-0.84458447) q[3];
sx q[3];
rz(-2.5648404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4429861) q[2];
sx q[2];
rz(-1.3977945) q[2];
sx q[2];
rz(-2.8473575) q[2];
rz(3.0596628) q[3];
sx q[3];
rz(-0.51920813) q[3];
sx q[3];
rz(3.1053655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0329523) q[0];
sx q[0];
rz(-2.2886798) q[0];
sx q[0];
rz(-0.0090573514) q[0];
rz(2.5065705) q[1];
sx q[1];
rz(-0.68990866) q[1];
sx q[1];
rz(-0.10805282) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2436284) q[0];
sx q[0];
rz(-1.122323) q[0];
sx q[0];
rz(1.9237706) q[0];
rz(1.7093541) q[2];
sx q[2];
rz(-1.1922622) q[2];
sx q[2];
rz(1.1059424) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6335771) q[1];
sx q[1];
rz(-1.0883691) q[1];
sx q[1];
rz(-2.0111994) q[1];
x q[2];
rz(-1.024527) q[3];
sx q[3];
rz(-1.3427991) q[3];
sx q[3];
rz(-2.8834164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.99284995) q[2];
sx q[2];
rz(-1.7603346) q[2];
sx q[2];
rz(-1.2711058) q[2];
rz(3.0631915) q[3];
sx q[3];
rz(-1.4784808) q[3];
sx q[3];
rz(1.1318077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8909797) q[0];
sx q[0];
rz(-1.5654634) q[0];
sx q[0];
rz(-0.65761956) q[0];
rz(1.744386) q[1];
sx q[1];
rz(-1.1746635) q[1];
sx q[1];
rz(0.89362842) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.052505715) q[0];
sx q[0];
rz(-2.1930709) q[0];
sx q[0];
rz(-1.1629348) q[0];
x q[1];
rz(-2.1003175) q[2];
sx q[2];
rz(-0.5404226) q[2];
sx q[2];
rz(1.4735917) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9295846) q[1];
sx q[1];
rz(-1.3699023) q[1];
sx q[1];
rz(-0.68540539) q[1];
rz(-pi) q[2];
rz(2.3818447) q[3];
sx q[3];
rz(-1.500251) q[3];
sx q[3];
rz(0.37978803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39067337) q[2];
sx q[2];
rz(-0.94736391) q[2];
sx q[2];
rz(0.52948362) q[2];
rz(0.47618619) q[3];
sx q[3];
rz(-1.809285) q[3];
sx q[3];
rz(2.7990394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7628409) q[0];
sx q[0];
rz(-1.6289926) q[0];
sx q[0];
rz(-1.09028) q[0];
rz(-0.11225637) q[1];
sx q[1];
rz(-1.1039762) q[1];
sx q[1];
rz(1.1539248) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40076462) q[0];
sx q[0];
rz(-2.4908998) q[0];
sx q[0];
rz(-3.0853737) q[0];
x q[1];
rz(2.5904028) q[2];
sx q[2];
rz(-1.4442208) q[2];
sx q[2];
rz(-1.7984496) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9824144) q[1];
sx q[1];
rz(-0.7632066) q[1];
sx q[1];
rz(3.0522507) q[1];
x q[2];
rz(-2.2418601) q[3];
sx q[3];
rz(-1.8814439) q[3];
sx q[3];
rz(1.336738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.551679) q[2];
sx q[2];
rz(-2.7605197) q[2];
sx q[2];
rz(2.7611458) q[2];
rz(1.1278661) q[3];
sx q[3];
rz(-1.7815855) q[3];
sx q[3];
rz(-1.5765223) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34148759) q[0];
sx q[0];
rz(-1.2416168) q[0];
sx q[0];
rz(-0.63968101) q[0];
rz(1.2387964) q[1];
sx q[1];
rz(-1.7265373) q[1];
sx q[1];
rz(-1.9715086) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1146678) q[0];
sx q[0];
rz(-1.2171193) q[0];
sx q[0];
rz(-0.35004079) q[0];
x q[1];
rz(0.77634546) q[2];
sx q[2];
rz(-1.3753969) q[2];
sx q[2];
rz(-2.5459144) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3635892) q[1];
sx q[1];
rz(-0.52007857) q[1];
sx q[1];
rz(-2.7705454) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3977259) q[3];
sx q[3];
rz(-0.67865463) q[3];
sx q[3];
rz(2.9242587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8273948) q[2];
sx q[2];
rz(-2.4908227) q[2];
sx q[2];
rz(0.38273746) q[2];
rz(-0.9283723) q[3];
sx q[3];
rz(-1.9675156) q[3];
sx q[3];
rz(-2.476957) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2255573) q[0];
sx q[0];
rz(-1.501843) q[0];
sx q[0];
rz(-2.8826707) q[0];
rz(-2.4312773) q[1];
sx q[1];
rz(-1.0537035) q[1];
sx q[1];
rz(2.6616667) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7326638) q[0];
sx q[0];
rz(-0.73584475) q[0];
sx q[0];
rz(-2.0793414) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9051022) q[2];
sx q[2];
rz(-0.91954008) q[2];
sx q[2];
rz(1.0054393) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.48206115) q[1];
sx q[1];
rz(-1.9226942) q[1];
sx q[1];
rz(2.5955276) q[1];
x q[2];
rz(-2.6142526) q[3];
sx q[3];
rz(-0.96685997) q[3];
sx q[3];
rz(-0.37249836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2991128) q[2];
sx q[2];
rz(-2.2128236) q[2];
sx q[2];
rz(1.8245565) q[2];
rz(-1.2420098) q[3];
sx q[3];
rz(-0.95359355) q[3];
sx q[3];
rz(0.73808134) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9823572) q[0];
sx q[0];
rz(-3.1095105) q[0];
sx q[0];
rz(1.6788917) q[0];
rz(-0.97933979) q[1];
sx q[1];
rz(-2.0420488) q[1];
sx q[1];
rz(2.2534823) q[1];
rz(-1.7469035) q[2];
sx q[2];
rz(-1.0955784) q[2];
sx q[2];
rz(-0.57217862) q[2];
rz(0.76673037) q[3];
sx q[3];
rz(-2.7646716) q[3];
sx q[3];
rz(-0.9203831) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
