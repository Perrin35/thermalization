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
rz(-0.55627745) q[0];
sx q[0];
rz(-0.039160691) q[0];
sx q[0];
rz(-0.66443366) q[0];
rz(2.9593664) q[1];
sx q[1];
rz(-1.6874977) q[1];
sx q[1];
rz(1.4092285) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8614757) q[0];
sx q[0];
rz(-1.404477) q[0];
sx q[0];
rz(-0.54027277) q[0];
rz(-1.9033236) q[2];
sx q[2];
rz(-1.9938139) q[2];
sx q[2];
rz(0.11655434) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.69610447) q[1];
sx q[1];
rz(-2.124212) q[1];
sx q[1];
rz(1.9457818) q[1];
rz(-2.9578311) q[3];
sx q[3];
rz(-1.2145044) q[3];
sx q[3];
rz(-0.7687062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0802143) q[2];
sx q[2];
rz(-2.5900216) q[2];
sx q[2];
rz(0.75458327) q[2];
rz(-1.4590229) q[3];
sx q[3];
rz(-2.2857234) q[3];
sx q[3];
rz(-1.4065019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.215312) q[0];
sx q[0];
rz(-0.54406852) q[0];
sx q[0];
rz(-1.1625483) q[0];
rz(-0.7112208) q[1];
sx q[1];
rz(-2.4237207) q[1];
sx q[1];
rz(-0.1444764) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9677749) q[0];
sx q[0];
rz(-1.4878232) q[0];
sx q[0];
rz(0.2936347) q[0];
rz(-pi) q[1];
rz(2.5108482) q[2];
sx q[2];
rz(-2.1877383) q[2];
sx q[2];
rz(-1.5396876) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.45119845) q[1];
sx q[1];
rz(-1.9393801) q[1];
sx q[1];
rz(-0.75779961) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9599592) q[3];
sx q[3];
rz(-0.70011052) q[3];
sx q[3];
rz(2.7896529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43313906) q[2];
sx q[2];
rz(-1.837156) q[2];
sx q[2];
rz(-0.30678314) q[2];
rz(2.3661546) q[3];
sx q[3];
rz(-0.38475761) q[3];
sx q[3];
rz(2.9966089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6590092) q[0];
sx q[0];
rz(-2.6982396) q[0];
sx q[0];
rz(-0.79202598) q[0];
rz(-1.5838985) q[1];
sx q[1];
rz(-1.3521103) q[1];
sx q[1];
rz(2.1477594) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4531012) q[0];
sx q[0];
rz(-1.9217886) q[0];
sx q[0];
rz(0.61601244) q[0];
rz(0.58980073) q[2];
sx q[2];
rz(-1.6720534) q[2];
sx q[2];
rz(-2.4112005) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6043338) q[1];
sx q[1];
rz(-0.52642614) q[1];
sx q[1];
rz(-0.98186214) q[1];
x q[2];
rz(-0.1900717) q[3];
sx q[3];
rz(-2.0160003) q[3];
sx q[3];
rz(-1.0560738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9516248) q[2];
sx q[2];
rz(-0.8364532) q[2];
sx q[2];
rz(0.17568976) q[2];
rz(-0.22819337) q[3];
sx q[3];
rz(-2.5722645) q[3];
sx q[3];
rz(0.5051676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2100385) q[0];
sx q[0];
rz(-2.2966972) q[0];
sx q[0];
rz(1.5616052) q[0];
rz(2.2700894) q[1];
sx q[1];
rz(-1.611064) q[1];
sx q[1];
rz(-0.16564381) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7270532) q[0];
sx q[0];
rz(-0.001984607) q[0];
sx q[0];
rz(-2.5059047) q[0];
rz(-pi) q[1];
x q[1];
rz(1.125096) q[2];
sx q[2];
rz(-1.5786849) q[2];
sx q[2];
rz(2.5119022) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6871275) q[1];
sx q[1];
rz(-1.4923959) q[1];
sx q[1];
rz(2.3509174) q[1];
rz(-2.5609162) q[3];
sx q[3];
rz(-1.3370974) q[3];
sx q[3];
rz(-0.068598821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.44276253) q[2];
sx q[2];
rz(-0.58219588) q[2];
sx q[2];
rz(-2.3637135) q[2];
rz(0.88464087) q[3];
sx q[3];
rz(-1.1529461) q[3];
sx q[3];
rz(2.1625429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1409461) q[0];
sx q[0];
rz(-2.1053173) q[0];
sx q[0];
rz(-2.272814) q[0];
rz(-2.2390305) q[1];
sx q[1];
rz(-1.9156009) q[1];
sx q[1];
rz(2.5696519) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8800655) q[0];
sx q[0];
rz(-2.2216317) q[0];
sx q[0];
rz(2.5291247) q[0];
x q[1];
rz(1.1210905) q[2];
sx q[2];
rz(-2.0313946) q[2];
sx q[2];
rz(3.1267303) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.16145591) q[1];
sx q[1];
rz(-1.2677578) q[1];
sx q[1];
rz(1.113446) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8945893) q[3];
sx q[3];
rz(-1.9387783) q[3];
sx q[3];
rz(-1.0447657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9945485) q[2];
sx q[2];
rz(-2.2037555) q[2];
sx q[2];
rz(-2.5082972) q[2];
rz(-0.58081943) q[3];
sx q[3];
rz(-1.7827026) q[3];
sx q[3];
rz(-0.66494989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0112515) q[0];
sx q[0];
rz(-2.7524152) q[0];
sx q[0];
rz(1.8036386) q[0];
rz(1.411571) q[1];
sx q[1];
rz(-0.4069702) q[1];
sx q[1];
rz(2.3445047) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10182285) q[0];
sx q[0];
rz(-1.2972141) q[0];
sx q[0];
rz(-2.9780529) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6930137) q[2];
sx q[2];
rz(-0.56899348) q[2];
sx q[2];
rz(-0.48556604) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6812207) q[1];
sx q[1];
rz(-2.8275194) q[1];
sx q[1];
rz(-1.0286147) q[1];
rz(-pi) q[2];
rz(1.0510526) q[3];
sx q[3];
rz(-1.4134348) q[3];
sx q[3];
rz(-1.5417772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.36593124) q[2];
sx q[2];
rz(-0.66880995) q[2];
sx q[2];
rz(-0.68525165) q[2];
rz(2.8819486) q[3];
sx q[3];
rz(-2.1681867) q[3];
sx q[3];
rz(-1.2118305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.892136) q[0];
sx q[0];
rz(-2.9849755) q[0];
sx q[0];
rz(0.53949612) q[0];
rz(2.9501713) q[1];
sx q[1];
rz(-1.6554662) q[1];
sx q[1];
rz(-0.44245455) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4726329) q[0];
sx q[0];
rz(-1.5711391) q[0];
sx q[0];
rz(1.5616722) q[0];
x q[1];
rz(-0.94687825) q[2];
sx q[2];
rz(-1.2420474) q[2];
sx q[2];
rz(-2.2524232) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4686582) q[1];
sx q[1];
rz(-0.71437049) q[1];
sx q[1];
rz(-2.0972392) q[1];
x q[2];
rz(-0.70480736) q[3];
sx q[3];
rz(-2.2205345) q[3];
sx q[3];
rz(-2.5832502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.29360867) q[2];
sx q[2];
rz(-0.30985761) q[2];
sx q[2];
rz(-0.7027182) q[2];
rz(-2.8637049) q[3];
sx q[3];
rz(-1.0669471) q[3];
sx q[3];
rz(-1.8028629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24669312) q[0];
sx q[0];
rz(-0.86240697) q[0];
sx q[0];
rz(-2.605751) q[0];
rz(2.4193173) q[1];
sx q[1];
rz(-1.4403789) q[1];
sx q[1];
rz(2.3586418) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7039258) q[0];
sx q[0];
rz(-1.9294318) q[0];
sx q[0];
rz(-3.0485247) q[0];
rz(-pi) q[1];
rz(-3.0083904) q[2];
sx q[2];
rz(-2.7410772) q[2];
sx q[2];
rz(-2.7513964) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6440575) q[1];
sx q[1];
rz(-0.78427197) q[1];
sx q[1];
rz(-2.3375744) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1184741) q[3];
sx q[3];
rz(-1.9503106) q[3];
sx q[3];
rz(1.4364157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7213514) q[2];
sx q[2];
rz(-2.6768117) q[2];
sx q[2];
rz(-2.6381524) q[2];
rz(0.33468801) q[3];
sx q[3];
rz(-1.3379593) q[3];
sx q[3];
rz(2.9065342) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47401416) q[0];
sx q[0];
rz(-0.36822167) q[0];
sx q[0];
rz(1.0598805) q[0];
rz(2.8267951) q[1];
sx q[1];
rz(-1.7960725) q[1];
sx q[1];
rz(1.559929) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81875694) q[0];
sx q[0];
rz(-2.4222932) q[0];
sx q[0];
rz(-2.8000205) q[0];
rz(-pi) q[1];
rz(1.8354206) q[2];
sx q[2];
rz(-1.6034063) q[2];
sx q[2];
rz(-0.6219686) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0421903) q[1];
sx q[1];
rz(-1.3563957) q[1];
sx q[1];
rz(2.2221788) q[1];
rz(-pi) q[2];
rz(-0.022059343) q[3];
sx q[3];
rz(-1.1556871) q[3];
sx q[3];
rz(2.4850415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9554837) q[2];
sx q[2];
rz(-1.9766221) q[2];
sx q[2];
rz(-1.0355518) q[2];
rz(-1.1459076) q[3];
sx q[3];
rz(-1.112554) q[3];
sx q[3];
rz(-0.15315332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.938852) q[0];
sx q[0];
rz(-0.11883141) q[0];
sx q[0];
rz(-2.3760997) q[0];
rz(1.0368404) q[1];
sx q[1];
rz(-1.8427269) q[1];
sx q[1];
rz(0.11229215) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2689026) q[0];
sx q[0];
rz(-1.2451474) q[0];
sx q[0];
rz(2.6652357) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61157558) q[2];
sx q[2];
rz(-2.2627352) q[2];
sx q[2];
rz(2.126978) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0794509) q[1];
sx q[1];
rz(-1.3704668) q[1];
sx q[1];
rz(-0.82605861) q[1];
rz(-pi) q[2];
x q[2];
rz(0.046979172) q[3];
sx q[3];
rz(-2.3571797) q[3];
sx q[3];
rz(1.7681153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3969193) q[2];
sx q[2];
rz(-1.8645218) q[2];
sx q[2];
rz(0.12167682) q[2];
rz(-2.7820898) q[3];
sx q[3];
rz(-2.4183141) q[3];
sx q[3];
rz(-3.1239037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59035463) q[0];
sx q[0];
rz(-1.4689162) q[0];
sx q[0];
rz(1.3015251) q[0];
rz(1.7474668) q[1];
sx q[1];
rz(-1.9448517) q[1];
sx q[1];
rz(1.3762884) q[1];
rz(2.5689498) q[2];
sx q[2];
rz(-2.4349458) q[2];
sx q[2];
rz(-1.4473421) q[2];
rz(-1.8471424) q[3];
sx q[3];
rz(-2.69097) q[3];
sx q[3];
rz(0.16926058) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
