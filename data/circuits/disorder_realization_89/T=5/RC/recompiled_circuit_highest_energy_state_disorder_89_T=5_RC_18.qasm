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
rz(1.8472449) q[0];
sx q[0];
rz(-1.6142774) q[0];
sx q[0];
rz(2.9543028) q[0];
rz(1.8118495) q[1];
sx q[1];
rz(-0.5464552) q[1];
sx q[1];
rz(1.7940914) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8719604) q[0];
sx q[0];
rz(-0.56817164) q[0];
sx q[0];
rz(-1.4251891) q[0];
rz(2.0894089) q[2];
sx q[2];
rz(-1.2402657) q[2];
sx q[2];
rz(1.7180819) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0653141) q[1];
sx q[1];
rz(-2.4597628) q[1];
sx q[1];
rz(1.5966459) q[1];
rz(-0.66708095) q[3];
sx q[3];
rz(-2.8192408) q[3];
sx q[3];
rz(-1.1548815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.61907855) q[2];
sx q[2];
rz(-1.9587025) q[2];
sx q[2];
rz(0.28156933) q[2];
rz(-1.36739) q[3];
sx q[3];
rz(-1.2308246) q[3];
sx q[3];
rz(-2.863133) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6047769) q[0];
sx q[0];
rz(-0.29412687) q[0];
sx q[0];
rz(-1.3500704) q[0];
rz(0.049909441) q[1];
sx q[1];
rz(-0.3388181) q[1];
sx q[1];
rz(1.2443589) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2771143) q[0];
sx q[0];
rz(-1.9196072) q[0];
sx q[0];
rz(-1.6326399) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2717541) q[2];
sx q[2];
rz(-0.56952945) q[2];
sx q[2];
rz(-2.6851252) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7283561) q[1];
sx q[1];
rz(-1.3750569) q[1];
sx q[1];
rz(1.2979862) q[1];
rz(-0.54173754) q[3];
sx q[3];
rz(-1.1835464) q[3];
sx q[3];
rz(-2.4335086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0656978) q[2];
sx q[2];
rz(-1.1827129) q[2];
sx q[2];
rz(0.069570216) q[2];
rz(-0.73404297) q[3];
sx q[3];
rz(-0.78248048) q[3];
sx q[3];
rz(-2.5231584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0760536) q[0];
sx q[0];
rz(-0.19936182) q[0];
sx q[0];
rz(-0.25654909) q[0];
rz(-1.8007295) q[1];
sx q[1];
rz(-2.1840265) q[1];
sx q[1];
rz(-2.0975013) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0437735) q[0];
sx q[0];
rz(-1.7390842) q[0];
sx q[0];
rz(1.8884482) q[0];
rz(0.759941) q[2];
sx q[2];
rz(-2.8360998) q[2];
sx q[2];
rz(-1.0585275) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.26790923) q[1];
sx q[1];
rz(-1.8716772) q[1];
sx q[1];
rz(1.8681225) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7875015) q[3];
sx q[3];
rz(-2.5575221) q[3];
sx q[3];
rz(1.3497373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9112245) q[2];
sx q[2];
rz(-2.0679097) q[2];
sx q[2];
rz(-0.71630859) q[2];
rz(-0.45799842) q[3];
sx q[3];
rz(-1.674998) q[3];
sx q[3];
rz(2.9108293) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3476747) q[0];
sx q[0];
rz(-2.1383998) q[0];
sx q[0];
rz(2.4804261) q[0];
rz(1.2541153) q[1];
sx q[1];
rz(-1.5690119) q[1];
sx q[1];
rz(1.5987781) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7044012) q[0];
sx q[0];
rz(-1.7462303) q[0];
sx q[0];
rz(-1.3153005) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44825596) q[2];
sx q[2];
rz(-1.0545306) q[2];
sx q[2];
rz(-0.22656952) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.99847163) q[1];
sx q[1];
rz(-2.4354092) q[1];
sx q[1];
rz(-0.097340214) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7518753) q[3];
sx q[3];
rz(-1.0709312) q[3];
sx q[3];
rz(3.0177649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8864112) q[2];
sx q[2];
rz(-2.7913783) q[2];
sx q[2];
rz(2.0070455) q[2];
rz(2.5875731) q[3];
sx q[3];
rz(-1.3628682) q[3];
sx q[3];
rz(-2.5567283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2968813) q[0];
sx q[0];
rz(-3.1394594) q[0];
sx q[0];
rz(-1.1312436) q[0];
rz(0.058111195) q[1];
sx q[1];
rz(-1.8986214) q[1];
sx q[1];
rz(1.4505454) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.03019985) q[0];
sx q[0];
rz(-0.95472017) q[0];
sx q[0];
rz(-2.6598823) q[0];
rz(-pi) q[1];
rz(1.8671473) q[2];
sx q[2];
rz(-0.83186281) q[2];
sx q[2];
rz(-0.30767627) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5202153) q[1];
sx q[1];
rz(-1.8978373) q[1];
sx q[1];
rz(-1.821063) q[1];
rz(-0.93673076) q[3];
sx q[3];
rz(-1.9104244) q[3];
sx q[3];
rz(-2.820211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6821373) q[2];
sx q[2];
rz(-2.5294999) q[2];
sx q[2];
rz(-1.3539782) q[2];
rz(2.5891384) q[3];
sx q[3];
rz(-0.8907291) q[3];
sx q[3];
rz(-0.49828211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8950997) q[0];
sx q[0];
rz(-2.1857078) q[0];
sx q[0];
rz(-1.0981052) q[0];
rz(1.8662628) q[1];
sx q[1];
rz(-1.9018491) q[1];
sx q[1];
rz(-1.0414418) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.784248) q[0];
sx q[0];
rz(-2.2929296) q[0];
sx q[0];
rz(-2.4289352) q[0];
rz(-2.8813754) q[2];
sx q[2];
rz(-0.12311664) q[2];
sx q[2];
rz(-2.3359063) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.66646092) q[1];
sx q[1];
rz(-2.0826) q[1];
sx q[1];
rz(-3.0252837) q[1];
x q[2];
rz(-1.249635) q[3];
sx q[3];
rz(-0.75157673) q[3];
sx q[3];
rz(2.1134714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2991403) q[2];
sx q[2];
rz(-1.1536529) q[2];
sx q[2];
rz(-0.057272043) q[2];
rz(0.20036571) q[3];
sx q[3];
rz(-0.55557957) q[3];
sx q[3];
rz(0.19181767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5361236) q[0];
sx q[0];
rz(-0.64490461) q[0];
sx q[0];
rz(2.7591925) q[0];
rz(2.1740225) q[1];
sx q[1];
rz(-0.41125527) q[1];
sx q[1];
rz(1.254902) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40036094) q[0];
sx q[0];
rz(-2.1810449) q[0];
sx q[0];
rz(-1.639354) q[0];
rz(-pi) q[1];
rz(-2.3284328) q[2];
sx q[2];
rz(-2.325468) q[2];
sx q[2];
rz(0.94674158) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.7130512) q[1];
sx q[1];
rz(-1.9425434) q[1];
sx q[1];
rz(2.0058238) q[1];
x q[2];
rz(-1.8462049) q[3];
sx q[3];
rz(-2.3304238) q[3];
sx q[3];
rz(1.0397183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3249698) q[2];
sx q[2];
rz(-2.7866456) q[2];
sx q[2];
rz(-0.85186446) q[2];
rz(-2.8999117) q[3];
sx q[3];
rz(-1.9342187) q[3];
sx q[3];
rz(0.91066256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8082751) q[0];
sx q[0];
rz(-0.18263826) q[0];
sx q[0];
rz(2.0727378) q[0];
rz(1.764027) q[1];
sx q[1];
rz(-2.865538) q[1];
sx q[1];
rz(-0.68881234) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23857164) q[0];
sx q[0];
rz(-1.8982783) q[0];
sx q[0];
rz(-0.37974289) q[0];
rz(0.19866069) q[2];
sx q[2];
rz(-2.5294494) q[2];
sx q[2];
rz(-2.8191301) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6530643) q[1];
sx q[1];
rz(-1.1875069) q[1];
sx q[1];
rz(-1.2252818) q[1];
rz(-0.28941713) q[3];
sx q[3];
rz(-2.5934813) q[3];
sx q[3];
rz(1.1374933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0869202) q[2];
sx q[2];
rz(-0.72267795) q[2];
sx q[2];
rz(1.5642536) q[2];
rz(-0.65586048) q[3];
sx q[3];
rz(-1.7003931) q[3];
sx q[3];
rz(-2.3724469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.61651388) q[0];
sx q[0];
rz(-1.0322796) q[0];
sx q[0];
rz(-2.9922564) q[0];
rz(-1.0859882) q[1];
sx q[1];
rz(-2.1859152) q[1];
sx q[1];
rz(0.48233262) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5457561) q[0];
sx q[0];
rz(-1.0991503) q[0];
sx q[0];
rz(-3.0455941) q[0];
x q[1];
rz(-1.3382089) q[2];
sx q[2];
rz(-1.5778189) q[2];
sx q[2];
rz(-2.1178031) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.73585549) q[1];
sx q[1];
rz(-2.2453826) q[1];
sx q[1];
rz(1.2184185) q[1];
rz(-0.41986044) q[3];
sx q[3];
rz(-1.4814988) q[3];
sx q[3];
rz(1.6557224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0444191) q[2];
sx q[2];
rz(-0.83690518) q[2];
sx q[2];
rz(-0.4591628) q[2];
rz(-3.0432213) q[3];
sx q[3];
rz(-1.8561074) q[3];
sx q[3];
rz(1.8927101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8654873) q[0];
sx q[0];
rz(-1.8806172) q[0];
sx q[0];
rz(-0.14044811) q[0];
rz(-1.2175995) q[1];
sx q[1];
rz(-2.0001037) q[1];
sx q[1];
rz(-1.5906895) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86027788) q[0];
sx q[0];
rz(-1.6522121) q[0];
sx q[0];
rz(2.415262) q[0];
x q[1];
rz(-1.7459007) q[2];
sx q[2];
rz(-1.4005454) q[2];
sx q[2];
rz(-2.424602) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5441011) q[1];
sx q[1];
rz(-0.95817503) q[1];
sx q[1];
rz(1.8922011) q[1];
rz(-pi) q[2];
rz(1.5545098) q[3];
sx q[3];
rz(-2.4706591) q[3];
sx q[3];
rz(0.85834225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6259674) q[2];
sx q[2];
rz(-2.1145861) q[2];
sx q[2];
rz(1.1930126) q[2];
rz(2.6595645) q[3];
sx q[3];
rz(-0.41976443) q[3];
sx q[3];
rz(1.0198063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-1.0057509) q[0];
sx q[0];
rz(-2.0477722) q[0];
sx q[0];
rz(2.3070106) q[0];
rz(-2.3781378) q[1];
sx q[1];
rz(-1.3044985) q[1];
sx q[1];
rz(0.34674092) q[1];
rz(3.1196567) q[2];
sx q[2];
rz(-1.4540945) q[2];
sx q[2];
rz(-3.1380063) q[2];
rz(1.9967358) q[3];
sx q[3];
rz(-1.0196232) q[3];
sx q[3];
rz(-2.41223) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
