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
rz(-1.3475013) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6997049) q[0];
sx q[0];
rz(-2.1322269) q[0];
sx q[0];
rz(-3.0492298) q[0];
x q[1];
rz(-2.7653469) q[2];
sx q[2];
rz(-2.0587789) q[2];
sx q[2];
rz(3.1057292) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.076278585) q[1];
sx q[1];
rz(-0.68182987) q[1];
sx q[1];
rz(-1.5449468) q[1];
x q[2];
rz(0.66708095) q[3];
sx q[3];
rz(-0.32235185) q[3];
sx q[3];
rz(-1.1548815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5225141) q[2];
sx q[2];
rz(-1.1828902) q[2];
sx q[2];
rz(-2.8600233) q[2];
rz(1.7742026) q[3];
sx q[3];
rz(-1.910768) q[3];
sx q[3];
rz(-0.27845964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6047769) q[0];
sx q[0];
rz(-2.8474658) q[0];
sx q[0];
rz(-1.7915223) q[0];
rz(0.049909441) q[1];
sx q[1];
rz(-0.3388181) q[1];
sx q[1];
rz(1.2443589) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2725213) q[0];
sx q[0];
rz(-1.5126813) q[0];
sx q[0];
rz(0.34942594) q[0];
rz(-pi) q[1];
rz(-1.0216818) q[2];
sx q[2];
rz(-1.7303338) q[2];
sx q[2];
rz(1.3683461) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4132365) q[1];
sx q[1];
rz(-1.3750569) q[1];
sx q[1];
rz(1.2979862) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66920993) q[3];
sx q[3];
rz(-2.4871181) q[3];
sx q[3];
rz(1.4231285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0758948) q[2];
sx q[2];
rz(-1.9588797) q[2];
sx q[2];
rz(-3.0720224) q[2];
rz(-2.4075497) q[3];
sx q[3];
rz(-0.78248048) q[3];
sx q[3];
rz(-0.61843425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.065539) q[0];
sx q[0];
rz(-2.9422308) q[0];
sx q[0];
rz(2.8850436) q[0];
rz(-1.8007295) q[1];
sx q[1];
rz(-2.1840265) q[1];
sx q[1];
rz(1.0440913) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0437735) q[0];
sx q[0];
rz(-1.7390842) q[0];
sx q[0];
rz(-1.2531444) q[0];
rz(2.3816517) q[2];
sx q[2];
rz(-2.8360998) q[2];
sx q[2];
rz(-2.0830652) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0716182) q[1];
sx q[1];
rz(-2.721792) q[1];
sx q[1];
rz(-2.384925) q[1];
x q[2];
rz(-0.99756156) q[3];
sx q[3];
rz(-1.4519534) q[3];
sx q[3];
rz(0.039418377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.23036817) q[2];
sx q[2];
rz(-1.0736829) q[2];
sx q[2];
rz(0.71630859) q[2];
rz(-2.6835942) q[3];
sx q[3];
rz(-1.4665946) q[3];
sx q[3];
rz(2.9108293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79391795) q[0];
sx q[0];
rz(-1.0031928) q[0];
sx q[0];
rz(-2.4804261) q[0];
rz(-1.2541153) q[1];
sx q[1];
rz(-1.5690119) q[1];
sx q[1];
rz(1.5428146) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9624299) q[0];
sx q[0];
rz(-1.822285) q[0];
sx q[0];
rz(2.9603987) q[0];
rz(-pi) q[1];
rz(-0.44825596) q[2];
sx q[2];
rz(-2.087062) q[2];
sx q[2];
rz(-2.9150231) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.49816541) q[1];
sx q[1];
rz(-1.6339059) q[1];
sx q[1];
rz(-2.4377512) q[1];
rz(-pi) q[2];
rz(-2.1786134) q[3];
sx q[3];
rz(-2.5180897) q[3];
sx q[3];
rz(-2.3092712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.25518146) q[2];
sx q[2];
rz(-0.35021439) q[2];
sx q[2];
rz(-1.1345471) q[2];
rz(0.55401951) q[3];
sx q[3];
rz(-1.3628682) q[3];
sx q[3];
rz(-0.58486432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2968813) q[0];
sx q[0];
rz(-0.0021332707) q[0];
sx q[0];
rz(2.010349) q[0];
rz(0.058111195) q[1];
sx q[1];
rz(-1.8986214) q[1];
sx q[1];
rz(-1.6910472) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7052225) q[0];
sx q[0];
rz(-2.3793325) q[0];
sx q[0];
rz(-0.99135474) q[0];
rz(-pi) q[1];
rz(-1.2744454) q[2];
sx q[2];
rz(-2.3097298) q[2];
sx q[2];
rz(-2.8339164) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.031353019) q[1];
sx q[1];
rz(-1.3340557) q[1];
sx q[1];
rz(-2.8048022) q[1];
rz(-pi) q[2];
rz(1.0330328) q[3];
sx q[3];
rz(-2.4335053) q[3];
sx q[3];
rz(1.6747812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6821373) q[2];
sx q[2];
rz(-0.61209279) q[2];
sx q[2];
rz(-1.7876145) q[2];
rz(-2.5891384) q[3];
sx q[3];
rz(-0.8907291) q[3];
sx q[3];
rz(0.49828211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.246493) q[0];
sx q[0];
rz(-0.95588481) q[0];
sx q[0];
rz(2.0434875) q[0];
rz(1.8662628) q[1];
sx q[1];
rz(-1.2397436) q[1];
sx q[1];
rz(1.0414418) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7324243) q[0];
sx q[0];
rz(-2.083626) q[0];
sx q[0];
rz(-0.70968117) q[0];
x q[1];
rz(-2.8813754) q[2];
sx q[2];
rz(-3.018476) q[2];
sx q[2];
rz(-0.80568635) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2944156) q[1];
sx q[1];
rz(-1.4694459) q[1];
sx q[1];
rz(1.0560929) q[1];
x q[2];
rz(-0.84534332) q[3];
sx q[3];
rz(-1.3535557) q[3];
sx q[3];
rz(-2.8373534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2991403) q[2];
sx q[2];
rz(-1.1536529) q[2];
sx q[2];
rz(0.057272043) q[2];
rz(0.20036571) q[3];
sx q[3];
rz(-0.55557957) q[3];
sx q[3];
rz(-2.949775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6054691) q[0];
sx q[0];
rz(-2.496688) q[0];
sx q[0];
rz(0.38240018) q[0];
rz(-0.96757013) q[1];
sx q[1];
rz(-2.7303374) q[1];
sx q[1];
rz(1.8866906) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28110958) q[0];
sx q[0];
rz(-2.5279928) q[0];
sx q[0];
rz(-3.0439418) q[0];
x q[1];
rz(0.81315984) q[2];
sx q[2];
rz(-2.325468) q[2];
sx q[2];
rz(-2.1948511) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6204024) q[1];
sx q[1];
rz(-2.5772021) q[1];
sx q[1];
rz(2.317313) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2953877) q[3];
sx q[3];
rz(-2.3304238) q[3];
sx q[3];
rz(1.0397183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3249698) q[2];
sx q[2];
rz(-0.35494706) q[2];
sx q[2];
rz(-2.2897282) q[2];
rz(0.24168092) q[3];
sx q[3];
rz(-1.9342187) q[3];
sx q[3];
rz(0.91066256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3333176) q[0];
sx q[0];
rz(-2.9589544) q[0];
sx q[0];
rz(1.0688548) q[0];
rz(-1.3775657) q[1];
sx q[1];
rz(-0.27605468) q[1];
sx q[1];
rz(0.68881234) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.903021) q[0];
sx q[0];
rz(-1.8982783) q[0];
sx q[0];
rz(-0.37974289) q[0];
rz(2.942932) q[2];
sx q[2];
rz(-2.5294494) q[2];
sx q[2];
rz(2.8191301) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.88676597) q[1];
sx q[1];
rz(-0.51023645) q[1];
sx q[1];
rz(-0.6986104) q[1];
x q[2];
rz(1.7433002) q[3];
sx q[3];
rz(-1.0478847) q[3];
sx q[3];
rz(2.3397818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0546725) q[2];
sx q[2];
rz(-0.72267795) q[2];
sx q[2];
rz(1.5773391) q[2];
rz(0.65586048) q[3];
sx q[3];
rz(-1.4411996) q[3];
sx q[3];
rz(-2.3724469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5250788) q[0];
sx q[0];
rz(-2.109313) q[0];
sx q[0];
rz(2.9922564) q[0];
rz(-2.0556045) q[1];
sx q[1];
rz(-2.1859152) q[1];
sx q[1];
rz(2.65926) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3369088) q[0];
sx q[0];
rz(-2.6609976) q[0];
sx q[0];
rz(-1.3850336) q[0];
rz(-pi) q[1];
rz(-3.1343757) q[2];
sx q[2];
rz(-1.8033779) q[2];
sx q[2];
rz(0.54534334) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0606976) q[1];
sx q[1];
rz(-1.8437064) q[1];
sx q[1];
rz(2.4358783) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7217322) q[3];
sx q[3];
rz(-1.6600939) q[3];
sx q[3];
rz(1.6557224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0444191) q[2];
sx q[2];
rz(-2.3046875) q[2];
sx q[2];
rz(-2.6824299) q[2];
rz(0.098371355) q[3];
sx q[3];
rz(-1.2854853) q[3];
sx q[3];
rz(1.2488825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8654873) q[0];
sx q[0];
rz(-1.2609755) q[0];
sx q[0];
rz(0.14044811) q[0];
rz(-1.2175995) q[1];
sx q[1];
rz(-1.1414889) q[1];
sx q[1];
rz(-1.5509031) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3397796) q[0];
sx q[0];
rz(-0.7300517) q[0];
sx q[0];
rz(-3.0193437) q[0];
x q[1];
rz(-0.17284278) q[2];
sx q[2];
rz(-1.3982492) q[2];
sx q[2];
rz(0.82383982) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1223476) q[1];
sx q[1];
rz(-2.4594895) q[1];
sx q[1];
rz(-2.7191799) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2416653) q[3];
sx q[3];
rz(-1.5809217) q[3];
sx q[3];
rz(2.4163818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6259674) q[2];
sx q[2];
rz(-2.1145861) q[2];
sx q[2];
rz(1.1930126) q[2];
rz(2.6595645) q[3];
sx q[3];
rz(-2.7218282) q[3];
sx q[3];
rz(-1.0198063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-1.7557549) q[2];
sx q[2];
rz(-3.0228563) q[2];
sx q[2];
rz(-2.9517662) q[2];
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
