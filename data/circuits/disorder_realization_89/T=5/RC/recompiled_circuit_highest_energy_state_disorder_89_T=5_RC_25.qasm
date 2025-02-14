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
rz(-1.3297431) q[1];
sx q[1];
rz(-2.5951374) q[1];
sx q[1];
rz(1.3475013) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8719604) q[0];
sx q[0];
rz(-2.573421) q[0];
sx q[0];
rz(-1.4251891) q[0];
rz(-0.37624575) q[2];
sx q[2];
rz(-2.0587789) q[2];
sx q[2];
rz(0.03586344) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.076278585) q[1];
sx q[1];
rz(-2.4597628) q[1];
sx q[1];
rz(1.5966459) q[1];
rz(-pi) q[2];
rz(0.25661664) q[3];
sx q[3];
rz(-1.3735176) q[3];
sx q[3];
rz(2.9159604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5225141) q[2];
sx q[2];
rz(-1.9587025) q[2];
sx q[2];
rz(2.8600233) q[2];
rz(1.7742026) q[3];
sx q[3];
rz(-1.2308246) q[3];
sx q[3];
rz(-2.863133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.6047769) q[0];
sx q[0];
rz(-2.8474658) q[0];
sx q[0];
rz(-1.3500704) q[0];
rz(-3.0916832) q[1];
sx q[1];
rz(-0.3388181) q[1];
sx q[1];
rz(-1.8972338) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8644783) q[0];
sx q[0];
rz(-1.9196072) q[0];
sx q[0];
rz(-1.6326399) q[0];
rz(-pi) q[1];
x q[1];
rz(0.18644615) q[2];
sx q[2];
rz(-1.0294339) q[2];
sx q[2];
rz(-0.10554927) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.10319732) q[1];
sx q[1];
rz(-1.8382676) q[1];
sx q[1];
rz(-2.9385376) q[1];
rz(-2.4723827) q[3];
sx q[3];
rz(-2.4871181) q[3];
sx q[3];
rz(-1.4231285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0656978) q[2];
sx q[2];
rz(-1.1827129) q[2];
sx q[2];
rz(3.0720224) q[2];
rz(2.4075497) q[3];
sx q[3];
rz(-0.78248048) q[3];
sx q[3];
rz(-2.5231584) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0760536) q[0];
sx q[0];
rz(-0.19936182) q[0];
sx q[0];
rz(2.8850436) q[0];
rz(1.8007295) q[1];
sx q[1];
rz(-2.1840265) q[1];
sx q[1];
rz(-1.0440913) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0437735) q[0];
sx q[0];
rz(-1.7390842) q[0];
sx q[0];
rz(1.2531444) q[0];
x q[1];
rz(0.22473904) q[2];
sx q[2];
rz(-1.3620952) q[2];
sx q[2];
rz(-1.2485742) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.26790923) q[1];
sx q[1];
rz(-1.8716772) q[1];
sx q[1];
rz(1.2734702) q[1];
rz(-pi) q[2];
rz(-0.14117853) q[3];
sx q[3];
rz(-1.0021082) q[3];
sx q[3];
rz(1.533826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.23036817) q[2];
sx q[2];
rz(-1.0736829) q[2];
sx q[2];
rz(2.4252841) q[2];
rz(0.45799842) q[3];
sx q[3];
rz(-1.674998) q[3];
sx q[3];
rz(-2.9108293) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79391795) q[0];
sx q[0];
rz(-2.1383998) q[0];
sx q[0];
rz(-2.4804261) q[0];
rz(1.8874774) q[1];
sx q[1];
rz(-1.5690119) q[1];
sx q[1];
rz(-1.5987781) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5970478) q[0];
sx q[0];
rz(-2.8327541) q[0];
sx q[0];
rz(-0.95914532) q[0];
rz(-0.44825596) q[2];
sx q[2];
rz(-1.0545306) q[2];
sx q[2];
rz(2.9150231) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.143121) q[1];
sx q[1];
rz(-2.4354092) q[1];
sx q[1];
rz(3.0442524) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1041252) q[3];
sx q[3];
rz(-1.9107483) q[3];
sx q[3];
rz(-1.8890017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.25518146) q[2];
sx q[2];
rz(-0.35021439) q[2];
sx q[2];
rz(-2.0070455) q[2];
rz(2.5875731) q[3];
sx q[3];
rz(-1.3628682) q[3];
sx q[3];
rz(-2.5567283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(1.2968813) q[0];
sx q[0];
rz(-0.0021332707) q[0];
sx q[0];
rz(-1.1312436) q[0];
rz(-3.0834815) q[1];
sx q[1];
rz(-1.2429712) q[1];
sx q[1];
rz(-1.4505454) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.03019985) q[0];
sx q[0];
rz(-0.95472017) q[0];
sx q[0];
rz(2.6598823) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3804316) q[2];
sx q[2];
rz(-1.3532172) q[2];
sx q[2];
rz(-2.0812931) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1102396) q[1];
sx q[1];
rz(-1.3340557) q[1];
sx q[1];
rz(-0.33679049) q[1];
x q[2];
rz(1.0330328) q[3];
sx q[3];
rz(-0.70808739) q[3];
sx q[3];
rz(-1.6747812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.45945534) q[2];
sx q[2];
rz(-0.61209279) q[2];
sx q[2];
rz(1.7876145) q[2];
rz(2.5891384) q[3];
sx q[3];
rz(-2.2508636) q[3];
sx q[3];
rz(-2.6433105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8950997) q[0];
sx q[0];
rz(-0.95588481) q[0];
sx q[0];
rz(2.0434875) q[0];
rz(-1.8662628) q[1];
sx q[1];
rz(-1.2397436) q[1];
sx q[1];
rz(-1.0414418) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35734468) q[0];
sx q[0];
rz(-0.84866306) q[0];
sx q[0];
rz(-0.71265744) q[0];
rz(-pi) q[1];
rz(0.11901151) q[2];
sx q[2];
rz(-1.5391943) q[2];
sx q[2];
rz(-1.0234444) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2944156) q[1];
sx q[1];
rz(-1.6721467) q[1];
sx q[1];
rz(-2.0854998) q[1];
x q[2];
rz(-1.249635) q[3];
sx q[3];
rz(-2.3900159) q[3];
sx q[3];
rz(1.0281212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2991403) q[2];
sx q[2];
rz(-1.1536529) q[2];
sx q[2];
rz(0.057272043) q[2];
rz(-0.20036571) q[3];
sx q[3];
rz(-2.5860131) q[3];
sx q[3];
rz(-2.949775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6054691) q[0];
sx q[0];
rz(-0.64490461) q[0];
sx q[0];
rz(-0.38240018) q[0];
rz(-2.1740225) q[1];
sx q[1];
rz(-2.7303374) q[1];
sx q[1];
rz(-1.8866906) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2097652) q[0];
sx q[0];
rz(-1.5146274) q[0];
sx q[0];
rz(2.5302391) q[0];
x q[1];
rz(-2.228565) q[2];
sx q[2];
rz(-1.0464729) q[2];
sx q[2];
rz(1.1990449) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5211902) q[1];
sx q[1];
rz(-2.5772021) q[1];
sx q[1];
rz(2.317313) q[1];
x q[2];
rz(-0.77882336) q[3];
sx q[3];
rz(-1.3723139) q[3];
sx q[3];
rz(-0.72328156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.81662285) q[2];
sx q[2];
rz(-2.7866456) q[2];
sx q[2];
rz(0.85186446) q[2];
rz(-2.8999117) q[3];
sx q[3];
rz(-1.207374) q[3];
sx q[3];
rz(2.2309301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3333176) q[0];
sx q[0];
rz(-0.18263826) q[0];
sx q[0];
rz(-1.0688548) q[0];
rz(-1.3775657) q[1];
sx q[1];
rz(-0.27605468) q[1];
sx q[1];
rz(0.68881234) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4599054) q[0];
sx q[0];
rz(-1.9294159) q[0];
sx q[0];
rz(1.2201409) q[0];
rz(-pi) q[1];
rz(0.19866069) q[2];
sx q[2];
rz(-2.5294494) q[2];
sx q[2];
rz(-2.8191301) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.090052) q[1];
sx q[1];
rz(-1.2512815) q[1];
sx q[1];
rz(0.40488147) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28941713) q[3];
sx q[3];
rz(-2.5934813) q[3];
sx q[3];
rz(-1.1374933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0546725) q[2];
sx q[2];
rz(-0.72267795) q[2];
sx q[2];
rz(1.5642536) q[2];
rz(0.65586048) q[3];
sx q[3];
rz(-1.4411996) q[3];
sx q[3];
rz(0.76914579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(2.5250788) q[0];
sx q[0];
rz(-1.0322796) q[0];
sx q[0];
rz(-0.14933625) q[0];
rz(-1.0859882) q[1];
sx q[1];
rz(-2.1859152) q[1];
sx q[1];
rz(0.48233262) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80468388) q[0];
sx q[0];
rz(-2.6609976) q[0];
sx q[0];
rz(-1.756559) q[0];
rz(-pi) q[1];
rz(-1.5403378) q[2];
sx q[2];
rz(-0.23269146) q[2];
sx q[2];
rz(-2.5649477) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.20374035) q[1];
sx q[1];
rz(-2.3934869) q[1];
sx q[1];
rz(-2.7341873) q[1];
x q[2];
rz(-0.21621426) q[3];
sx q[3];
rz(-0.42869854) q[3];
sx q[3];
rz(-0.11224953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0444191) q[2];
sx q[2];
rz(-2.3046875) q[2];
sx q[2];
rz(0.4591628) q[2];
rz(-3.0432213) q[3];
sx q[3];
rz(-1.2854853) q[3];
sx q[3];
rz(1.2488825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2761053) q[0];
sx q[0];
rz(-1.2609755) q[0];
sx q[0];
rz(3.0011445) q[0];
rz(1.9239931) q[1];
sx q[1];
rz(-1.1414889) q[1];
sx q[1];
rz(-1.5509031) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5031932) q[0];
sx q[0];
rz(-0.84740438) q[0];
sx q[0];
rz(1.6795078) q[0];
x q[1];
rz(-2.9687499) q[2];
sx q[2];
rz(-1.3982492) q[2];
sx q[2];
rz(-0.82383982) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.019245) q[1];
sx q[1];
rz(-2.4594895) q[1];
sx q[1];
rz(-0.42241272) q[1];
x q[2];
rz(3.1286661) q[3];
sx q[3];
rz(-0.89996808) q[3];
sx q[3];
rz(-2.304043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6259674) q[2];
sx q[2];
rz(-1.0270065) q[2];
sx q[2];
rz(1.1930126) q[2];
rz(-0.48202816) q[3];
sx q[3];
rz(-0.41976443) q[3];
sx q[3];
rz(1.0198063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0057509) q[0];
sx q[0];
rz(-1.0938205) q[0];
sx q[0];
rz(-0.83458207) q[0];
rz(2.3781378) q[1];
sx q[1];
rz(-1.8370942) q[1];
sx q[1];
rz(-2.7948517) q[1];
rz(1.7557549) q[2];
sx q[2];
rz(-0.11873636) q[2];
sx q[2];
rz(0.18982646) q[2];
rz(1.1448568) q[3];
sx q[3];
rz(-2.1219694) q[3];
sx q[3];
rz(0.72936264) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
