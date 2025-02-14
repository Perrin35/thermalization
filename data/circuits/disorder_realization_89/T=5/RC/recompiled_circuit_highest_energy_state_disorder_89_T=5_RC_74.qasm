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
rz(1.5273153) q[0];
sx q[0];
rz(9.6120678) q[0];
rz(1.8118495) q[1];
sx q[1];
rz(-0.5464552) q[1];
sx q[1];
rz(-1.3475013) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8719604) q[0];
sx q[0];
rz(-0.56817164) q[0];
sx q[0];
rz(-1.7164036) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1762835) q[2];
sx q[2];
rz(-0.60676736) q[2];
sx q[2];
rz(-2.4772522) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.076278585) q[1];
sx q[1];
rz(-0.68182987) q[1];
sx q[1];
rz(-1.5449468) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7745733) q[3];
sx q[3];
rz(-1.319266) q[3];
sx q[3];
rz(-1.8478145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5225141) q[2];
sx q[2];
rz(-1.9587025) q[2];
sx q[2];
rz(-0.28156933) q[2];
rz(1.7742026) q[3];
sx q[3];
rz(-1.910768) q[3];
sx q[3];
rz(-0.27845964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53681579) q[0];
sx q[0];
rz(-2.8474658) q[0];
sx q[0];
rz(-1.3500704) q[0];
rz(-3.0916832) q[1];
sx q[1];
rz(-2.8027746) q[1];
sx q[1];
rz(1.8972338) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6852463) q[0];
sx q[0];
rz(-2.7875627) q[0];
sx q[0];
rz(-0.1683321) q[0];
x q[1];
rz(-1.0216818) q[2];
sx q[2];
rz(-1.4112588) q[2];
sx q[2];
rz(1.7732466) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7283561) q[1];
sx q[1];
rz(-1.3750569) q[1];
sx q[1];
rz(1.8436064) q[1];
rz(-pi) q[2];
rz(-1.1265321) q[3];
sx q[3];
rz(-1.0730529) q[3];
sx q[3];
rz(2.5023482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0758948) q[2];
sx q[2];
rz(-1.9588797) q[2];
sx q[2];
rz(-0.069570216) q[2];
rz(-0.73404297) q[3];
sx q[3];
rz(-2.3591122) q[3];
sx q[3];
rz(-0.61843425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.065539) q[0];
sx q[0];
rz(-2.9422308) q[0];
sx q[0];
rz(-2.8850436) q[0];
rz(-1.8007295) q[1];
sx q[1];
rz(-2.1840265) q[1];
sx q[1];
rz(1.0440913) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.097819177) q[0];
sx q[0];
rz(-1.4025084) q[0];
sx q[0];
rz(-1.2531444) q[0];
x q[1];
rz(-1.7847204) q[2];
sx q[2];
rz(-1.3510149) q[2];
sx q[2];
rz(0.27489382) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0699745) q[1];
sx q[1];
rz(-0.41980068) q[1];
sx q[1];
rz(-2.384925) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1440311) q[3];
sx q[3];
rz(-1.6896393) q[3];
sx q[3];
rz(-3.1021743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9112245) q[2];
sx q[2];
rz(-1.0736829) q[2];
sx q[2];
rz(0.71630859) q[2];
rz(0.45799842) q[3];
sx q[3];
rz(-1.4665946) q[3];
sx q[3];
rz(2.9108293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-0.79391795) q[0];
sx q[0];
rz(-1.0031928) q[0];
sx q[0];
rz(-2.4804261) q[0];
rz(1.8874774) q[1];
sx q[1];
rz(-1.5725807) q[1];
sx q[1];
rz(1.5987781) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4371915) q[0];
sx q[0];
rz(-1.7462303) q[0];
sx q[0];
rz(1.3153005) q[0];
rz(-pi) q[1];
x q[1];
rz(0.44825596) q[2];
sx q[2];
rz(-1.0545306) q[2];
sx q[2];
rz(-2.9150231) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6434272) q[1];
sx q[1];
rz(-1.6339059) q[1];
sx q[1];
rz(0.7038415) q[1];
rz(1.0374674) q[3];
sx q[3];
rz(-1.9107483) q[3];
sx q[3];
rz(1.8890017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.25518146) q[2];
sx q[2];
rz(-2.7913783) q[2];
sx q[2];
rz(-2.0070455) q[2];
rz(0.55401951) q[3];
sx q[3];
rz(-1.7787245) q[3];
sx q[3];
rz(-2.5567283) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8447113) q[0];
sx q[0];
rz(-3.1394594) q[0];
sx q[0];
rz(-1.1312436) q[0];
rz(-0.058111195) q[1];
sx q[1];
rz(-1.8986214) q[1];
sx q[1];
rz(-1.4505454) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.03019985) q[0];
sx q[0];
rz(-2.1868725) q[0];
sx q[0];
rz(-2.6598823) q[0];
rz(2.8314231) q[2];
sx q[2];
rz(-2.355994) q[2];
sx q[2];
rz(-2.4082662) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5202153) q[1];
sx q[1];
rz(-1.8978373) q[1];
sx q[1];
rz(-1.821063) q[1];
rz(-pi) q[2];
rz(-2.1085598) q[3];
sx q[3];
rz(-0.70808739) q[3];
sx q[3];
rz(1.4668114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6821373) q[2];
sx q[2];
rz(-2.5294999) q[2];
sx q[2];
rz(-1.7876145) q[2];
rz(-0.55245429) q[3];
sx q[3];
rz(-2.2508636) q[3];
sx q[3];
rz(-2.6433105) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8950997) q[0];
sx q[0];
rz(-2.1857078) q[0];
sx q[0];
rz(-2.0434875) q[0];
rz(-1.2753298) q[1];
sx q[1];
rz(-1.9018491) q[1];
sx q[1];
rz(2.1001508) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.784248) q[0];
sx q[0];
rz(-0.84866306) q[0];
sx q[0];
rz(0.71265744) q[0];
x q[1];
rz(1.5389693) q[2];
sx q[2];
rz(-1.4518445) q[2];
sx q[2];
rz(-2.5980189) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4751317) q[1];
sx q[1];
rz(-1.0589927) q[1];
sx q[1];
rz(0.11630897) q[1];
x q[2];
rz(-2.8547229) q[3];
sx q[3];
rz(-2.2755945) q[3];
sx q[3];
rz(-1.4554086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8424524) q[2];
sx q[2];
rz(-1.9879397) q[2];
sx q[2];
rz(-3.0843206) q[2];
rz(-2.9412269) q[3];
sx q[3];
rz(-0.55557957) q[3];
sx q[3];
rz(0.19181767) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6054691) q[0];
sx q[0];
rz(-0.64490461) q[0];
sx q[0];
rz(2.7591925) q[0];
rz(2.1740225) q[1];
sx q[1];
rz(-0.41125527) q[1];
sx q[1];
rz(1.254902) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7412317) q[0];
sx q[0];
rz(-0.96054777) q[0];
sx q[0];
rz(1.639354) q[0];
rz(-0.91302769) q[2];
sx q[2];
rz(-2.0951197) q[2];
sx q[2];
rz(1.1990449) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6204024) q[1];
sx q[1];
rz(-0.56439059) q[1];
sx q[1];
rz(-0.82427967) q[1];
x q[2];
rz(0.27886919) q[3];
sx q[3];
rz(-0.79852102) q[3];
sx q[3];
rz(2.4912733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3249698) q[2];
sx q[2];
rz(-2.7866456) q[2];
sx q[2];
rz(2.2897282) q[2];
rz(-0.24168092) q[3];
sx q[3];
rz(-1.9342187) q[3];
sx q[3];
rz(-0.91066256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3333176) q[0];
sx q[0];
rz(-0.18263826) q[0];
sx q[0];
rz(1.0688548) q[0];
rz(1.3775657) q[1];
sx q[1];
rz(-0.27605468) q[1];
sx q[1];
rz(2.4527803) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23857164) q[0];
sx q[0];
rz(-1.2433143) q[0];
sx q[0];
rz(2.7618498) q[0];
rz(-pi) q[1];
rz(-0.19866069) q[2];
sx q[2];
rz(-0.61214329) q[2];
sx q[2];
rz(-2.8191301) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.88676597) q[1];
sx q[1];
rz(-0.51023645) q[1];
sx q[1];
rz(-0.6986104) q[1];
rz(-pi) q[2];
rz(0.28941713) q[3];
sx q[3];
rz(-2.5934813) q[3];
sx q[3];
rz(2.0040994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0869202) q[2];
sx q[2];
rz(-0.72267795) q[2];
sx q[2];
rz(-1.5642536) q[2];
rz(-2.4857322) q[3];
sx q[3];
rz(-1.4411996) q[3];
sx q[3];
rz(0.76914579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61651388) q[0];
sx q[0];
rz(-1.0322796) q[0];
sx q[0];
rz(-2.9922564) q[0];
rz(-2.0556045) q[1];
sx q[1];
rz(-2.1859152) q[1];
sx q[1];
rz(-0.48233262) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80468388) q[0];
sx q[0];
rz(-2.6609976) q[0];
sx q[0];
rz(-1.756559) q[0];
rz(-pi) q[1];
rz(1.5403378) q[2];
sx q[2];
rz(-2.9089012) q[2];
sx q[2];
rz(-2.5649477) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4057372) q[1];
sx q[1];
rz(-2.2453826) q[1];
sx q[1];
rz(1.9231741) q[1];
rz(1.4730569) q[3];
sx q[3];
rz(-1.1527138) q[3];
sx q[3];
rz(-0.12471499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.09717354) q[2];
sx q[2];
rz(-0.83690518) q[2];
sx q[2];
rz(-0.4591628) q[2];
rz(3.0432213) q[3];
sx q[3];
rz(-1.8561074) q[3];
sx q[3];
rz(1.2488825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2761053) q[0];
sx q[0];
rz(-1.2609755) q[0];
sx q[0];
rz(-3.0011445) q[0];
rz(1.9239931) q[1];
sx q[1];
rz(-2.0001037) q[1];
sx q[1];
rz(-1.5906895) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2813148) q[0];
sx q[0];
rz(-1.4893805) q[0];
sx q[0];
rz(-0.72633065) q[0];
rz(0.17284278) q[2];
sx q[2];
rz(-1.3982492) q[2];
sx q[2];
rz(-0.82383982) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.925736) q[1];
sx q[1];
rz(-1.832215) q[1];
sx q[1];
rz(2.5040313) q[1];
rz(-pi) q[2];
rz(0.89992739) q[3];
sx q[3];
rz(-1.560671) q[3];
sx q[3];
rz(-2.4163818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5156252) q[2];
sx q[2];
rz(-1.0270065) q[2];
sx q[2];
rz(-1.94858) q[2];
rz(0.48202816) q[3];
sx q[3];
rz(-0.41976443) q[3];
sx q[3];
rz(2.1217864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-1.3858377) q[2];
sx q[2];
rz(-0.11873636) q[2];
sx q[2];
rz(0.18982646) q[2];
rz(-1.9967358) q[3];
sx q[3];
rz(-2.1219694) q[3];
sx q[3];
rz(0.72936264) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
