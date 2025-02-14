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
rz(2.5853152) q[0];
sx q[0];
rz(-3.102432) q[0];
sx q[0];
rz(-2.477159) q[0];
rz(-0.18222624) q[1];
sx q[1];
rz(-1.4540949) q[1];
sx q[1];
rz(1.7323642) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0213703) q[0];
sx q[0];
rz(-2.5787368) q[0];
sx q[0];
rz(-2.8261306) q[0];
x q[1];
rz(1.9033236) q[2];
sx q[2];
rz(-1.1477787) q[2];
sx q[2];
rz(-3.0250383) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.69610447) q[1];
sx q[1];
rz(-2.124212) q[1];
sx q[1];
rz(-1.1958108) q[1];
x q[2];
rz(-1.9326747) q[3];
sx q[3];
rz(-1.7428977) q[3];
sx q[3];
rz(-0.73735305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0802143) q[2];
sx q[2];
rz(-2.5900216) q[2];
sx q[2];
rz(0.75458327) q[2];
rz(-1.6825698) q[3];
sx q[3];
rz(-0.85586923) q[3];
sx q[3];
rz(-1.4065019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9262806) q[0];
sx q[0];
rz(-2.5975241) q[0];
sx q[0];
rz(-1.1625483) q[0];
rz(2.4303719) q[1];
sx q[1];
rz(-0.7178719) q[1];
sx q[1];
rz(0.1444764) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9677749) q[0];
sx q[0];
rz(-1.4878232) q[0];
sx q[0];
rz(-0.2936347) q[0];
x q[1];
rz(-2.2644193) q[2];
sx q[2];
rz(-0.85169221) q[2];
sx q[2];
rz(-2.4404877) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6934168) q[1];
sx q[1];
rz(-0.87478335) q[1];
sx q[1];
rz(-2.059518) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8322093) q[3];
sx q[3];
rz(-0.9321292) q[3];
sx q[3];
rz(2.2974599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.43313906) q[2];
sx q[2];
rz(-1.837156) q[2];
sx q[2];
rz(-0.30678314) q[2];
rz(2.3661546) q[3];
sx q[3];
rz(-2.756835) q[3];
sx q[3];
rz(0.14498372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6590092) q[0];
sx q[0];
rz(-0.44335303) q[0];
sx q[0];
rz(-0.79202598) q[0];
rz(1.5576942) q[1];
sx q[1];
rz(-1.7894824) q[1];
sx q[1];
rz(-2.1477594) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8071459) q[0];
sx q[0];
rz(-2.4440571) q[0];
sx q[0];
rz(2.5767482) q[0];
x q[1];
rz(2.9609072) q[2];
sx q[2];
rz(-0.59741086) q[2];
sx q[2];
rz(2.1512845) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0210304) q[1];
sx q[1];
rz(-1.1397727) q[1];
sx q[1];
rz(2.8293508) q[1];
rz(-pi) q[2];
rz(-2.023104) q[3];
sx q[3];
rz(-1.399446) q[3];
sx q[3];
rz(-2.7095344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1899679) q[2];
sx q[2];
rz(-2.3051395) q[2];
sx q[2];
rz(0.17568976) q[2];
rz(2.9133993) q[3];
sx q[3];
rz(-0.56932813) q[3];
sx q[3];
rz(2.6364251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2100385) q[0];
sx q[0];
rz(-2.2966972) q[0];
sx q[0];
rz(1.5799874) q[0];
rz(-0.87150323) q[1];
sx q[1];
rz(-1.611064) q[1];
sx q[1];
rz(-0.16564381) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3627421) q[0];
sx q[0];
rz(-1.5691994) q[0];
sx q[0];
rz(-1.569618) q[0];
rz(-pi) q[1];
rz(-1.5890938) q[2];
sx q[2];
rz(-0.4457655) q[2];
sx q[2];
rz(2.183977) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1936744) q[1];
sx q[1];
rz(-2.3478824) q[1];
sx q[1];
rz(0.11007424) q[1];
rz(0.40940855) q[3];
sx q[3];
rz(-0.62088517) q[3];
sx q[3];
rz(-1.3001022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.44276253) q[2];
sx q[2];
rz(-2.5593968) q[2];
sx q[2];
rz(0.77787918) q[2];
rz(-2.2569518) q[3];
sx q[3];
rz(-1.1529461) q[3];
sx q[3];
rz(2.1625429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1409461) q[0];
sx q[0];
rz(-1.0362754) q[0];
sx q[0];
rz(2.272814) q[0];
rz(-2.2390305) q[1];
sx q[1];
rz(-1.2259918) q[1];
sx q[1];
rz(0.57194078) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2615271) q[0];
sx q[0];
rz(-0.91996096) q[0];
sx q[0];
rz(-2.5291247) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1210905) q[2];
sx q[2];
rz(-1.110198) q[2];
sx q[2];
rz(3.1267303) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.16145591) q[1];
sx q[1];
rz(-1.8738349) q[1];
sx q[1];
rz(1.113446) q[1];
rz(-pi) q[2];
rz(2.8945893) q[3];
sx q[3];
rz(-1.2028143) q[3];
sx q[3];
rz(-1.0447657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9945485) q[2];
sx q[2];
rz(-2.2037555) q[2];
sx q[2];
rz(2.5082972) q[2];
rz(-0.58081943) q[3];
sx q[3];
rz(-1.7827026) q[3];
sx q[3];
rz(2.4766428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-3.0112515) q[0];
sx q[0];
rz(-2.7524152) q[0];
sx q[0];
rz(1.8036386) q[0];
rz(1.7300216) q[1];
sx q[1];
rz(-0.4069702) q[1];
sx q[1];
rz(0.79708797) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4244193) q[0];
sx q[0];
rz(-1.4133905) q[0];
sx q[0];
rz(-1.2936997) q[0];
rz(-pi) q[1];
x q[1];
rz(1.448579) q[2];
sx q[2];
rz(-0.56899348) q[2];
sx q[2];
rz(0.48556604) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.89576777) q[1];
sx q[1];
rz(-1.8386158) q[1];
sx q[1];
rz(-0.16606776) q[1];
rz(-pi) q[2];
rz(-1.0510526) q[3];
sx q[3];
rz(-1.7281579) q[3];
sx q[3];
rz(1.5998154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7756614) q[2];
sx q[2];
rz(-2.4727827) q[2];
sx q[2];
rz(2.456341) q[2];
rz(-0.25964409) q[3];
sx q[3];
rz(-2.1681867) q[3];
sx q[3];
rz(1.9297622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
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
rz(2.6991381) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.060613077) q[0];
sx q[0];
rz(-0.0091305841) q[0];
sx q[0];
rz(-1.5332444) q[0];
rz(-pi) q[1];
rz(0.94687825) q[2];
sx q[2];
rz(-1.8995452) q[2];
sx q[2];
rz(-2.2524232) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.67293447) q[1];
sx q[1];
rz(-0.71437049) q[1];
sx q[1];
rz(-1.0443535) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4367853) q[3];
sx q[3];
rz(-0.92105812) q[3];
sx q[3];
rz(0.55834246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.847984) q[2];
sx q[2];
rz(-2.831735) q[2];
sx q[2];
rz(-2.4388745) q[2];
rz(0.27788776) q[3];
sx q[3];
rz(-2.0746456) q[3];
sx q[3];
rz(1.8028629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24669312) q[0];
sx q[0];
rz(-2.2791857) q[0];
sx q[0];
rz(2.605751) q[0];
rz(0.72227532) q[1];
sx q[1];
rz(-1.4403789) q[1];
sx q[1];
rz(0.78295082) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43766681) q[0];
sx q[0];
rz(-1.2121608) q[0];
sx q[0];
rz(-0.093067972) q[0];
x q[1];
rz(3.0083904) q[2];
sx q[2];
rz(-2.7410772) q[2];
sx q[2];
rz(2.7513964) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4346501) q[1];
sx q[1];
rz(-1.0371814) q[1];
sx q[1];
rz(-0.60551079) q[1];
rz(2.1184741) q[3];
sx q[3];
rz(-1.9503106) q[3];
sx q[3];
rz(-1.4364157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7213514) q[2];
sx q[2];
rz(-2.6768117) q[2];
sx q[2];
rz(-0.50344023) q[2];
rz(-0.33468801) q[3];
sx q[3];
rz(-1.8036333) q[3];
sx q[3];
rz(2.9065342) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6675785) q[0];
sx q[0];
rz(-0.36822167) q[0];
sx q[0];
rz(1.0598805) q[0];
rz(-2.8267951) q[1];
sx q[1];
rz(-1.7960725) q[1];
sx q[1];
rz(1.5816636) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8813635) q[0];
sx q[0];
rz(-0.90103982) q[0];
sx q[0];
rz(-1.8561646) q[0];
rz(-pi) q[1];
rz(-1.8354206) q[2];
sx q[2];
rz(-1.6034063) q[2];
sx q[2];
rz(0.6219686) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4521802) q[1];
sx q[1];
rz(-2.204837) q[1];
sx q[1];
rz(0.26726064) q[1];
x q[2];
rz(-1.1555973) q[3];
sx q[3];
rz(-1.590982) q[3];
sx q[3];
rz(2.2184499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9554837) q[2];
sx q[2];
rz(-1.1649705) q[2];
sx q[2];
rz(2.1060409) q[2];
rz(-1.995685) q[3];
sx q[3];
rz(-1.112554) q[3];
sx q[3];
rz(0.15315332) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.938852) q[0];
sx q[0];
rz(-3.0227612) q[0];
sx q[0];
rz(0.76549292) q[0];
rz(-1.0368404) q[1];
sx q[1];
rz(-1.8427269) q[1];
sx q[1];
rz(3.0293005) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.87269) q[0];
sx q[0];
rz(-1.2451474) q[0];
sx q[0];
rz(-2.6652357) q[0];
rz(-pi) q[1];
rz(-2.3621777) q[2];
sx q[2];
rz(-2.0287435) q[2];
sx q[2];
rz(-0.13546695) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0621418) q[1];
sx q[1];
rz(-1.3704668) q[1];
sx q[1];
rz(0.82605861) q[1];
rz(0.046979172) q[3];
sx q[3];
rz(-2.3571797) q[3];
sx q[3];
rz(1.7681153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3969193) q[2];
sx q[2];
rz(-1.8645218) q[2];
sx q[2];
rz(0.12167682) q[2];
rz(-2.7820898) q[3];
sx q[3];
rz(-2.4183141) q[3];
sx q[3];
rz(0.017688964) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59035463) q[0];
sx q[0];
rz(-1.4689162) q[0];
sx q[0];
rz(1.3015251) q[0];
rz(1.3941258) q[1];
sx q[1];
rz(-1.196741) q[1];
sx q[1];
rz(-1.7653042) q[1];
rz(-2.5192026) q[2];
sx q[2];
rz(-1.211282) q[2];
sx q[2];
rz(-0.33242339) q[2];
rz(0.13124851) q[3];
sx q[3];
rz(-2.003142) q[3];
sx q[3];
rz(3.005645) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
