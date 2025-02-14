OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5366323) q[0];
sx q[0];
rz(-0.16083117) q[0];
sx q[0];
rz(2.7503126) q[0];
rz(3.0985576) q[1];
sx q[1];
rz(-1.6208836) q[1];
sx q[1];
rz(0.33837858) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94253263) q[0];
sx q[0];
rz(-1.8359658) q[0];
sx q[0];
rz(1.8607748) q[0];
x q[1];
rz(3.0909371) q[2];
sx q[2];
rz(-1.7850375) q[2];
sx q[2];
rz(3.0681477) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4520054) q[1];
sx q[1];
rz(-0.70065672) q[1];
sx q[1];
rz(-1.868471) q[1];
rz(-pi) q[2];
x q[2];
rz(1.67946) q[3];
sx q[3];
rz(-1.9851284) q[3];
sx q[3];
rz(-1.3760096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7732064) q[2];
sx q[2];
rz(-2.1619004) q[2];
sx q[2];
rz(2.0835908) q[2];
rz(0.10231415) q[3];
sx q[3];
rz(-1.5240069) q[3];
sx q[3];
rz(-1.4003096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31829396) q[0];
sx q[0];
rz(-0.75495356) q[0];
sx q[0];
rz(-2.9276983) q[0];
rz(1.8077883) q[1];
sx q[1];
rz(-2.6413481) q[1];
sx q[1];
rz(-2.852829) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64015111) q[0];
sx q[0];
rz(-2.1955865) q[0];
sx q[0];
rz(0.21296091) q[0];
rz(-pi) q[1];
rz(2.2877588) q[2];
sx q[2];
rz(-0.87426155) q[2];
sx q[2];
rz(2.8643212) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.60464478) q[1];
sx q[1];
rz(-1.5411475) q[1];
sx q[1];
rz(1.7916337) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2330194) q[3];
sx q[3];
rz(-1.1572517) q[3];
sx q[3];
rz(3.0440083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5220962) q[2];
sx q[2];
rz(-2.3213826) q[2];
sx q[2];
rz(-1.8572469) q[2];
rz(1.8340825) q[3];
sx q[3];
rz(-1.3176354) q[3];
sx q[3];
rz(-2.4431084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.4035325) q[0];
sx q[0];
rz(-2.0750676) q[0];
sx q[0];
rz(2.2962978) q[0];
rz(-2.2360133) q[1];
sx q[1];
rz(-2.5366492) q[1];
sx q[1];
rz(-1.9639429) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2098006) q[0];
sx q[0];
rz(-1.5508979) q[0];
sx q[0];
rz(-0.09505247) q[0];
rz(-pi) q[1];
rz(-0.81434701) q[2];
sx q[2];
rz(-0.6993412) q[2];
sx q[2];
rz(-1.7244848) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2237071) q[1];
sx q[1];
rz(-1.0038923) q[1];
sx q[1];
rz(-1.8170243) q[1];
rz(-pi) q[2];
rz(2.2159141) q[3];
sx q[3];
rz(-1.698602) q[3];
sx q[3];
rz(-2.5296581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5251069) q[2];
sx q[2];
rz(-0.944204) q[2];
sx q[2];
rz(3.0109829) q[2];
rz(-1.3579926) q[3];
sx q[3];
rz(-1.0087174) q[3];
sx q[3];
rz(1.85359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4865049) q[0];
sx q[0];
rz(-2.8772652) q[0];
sx q[0];
rz(0.41410145) q[0];
rz(-2.2762903) q[1];
sx q[1];
rz(-1.9046141) q[1];
sx q[1];
rz(-0.6032595) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8651486) q[0];
sx q[0];
rz(-1.28946) q[0];
sx q[0];
rz(1.9758609) q[0];
x q[1];
rz(-1.6418936) q[2];
sx q[2];
rz(-1.9219766) q[2];
sx q[2];
rz(-0.49654135) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.328426) q[1];
sx q[1];
rz(-0.35086497) q[1];
sx q[1];
rz(-0.12520949) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2141718) q[3];
sx q[3];
rz(-2.8093336) q[3];
sx q[3];
rz(1.3300524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.52411756) q[2];
sx q[2];
rz(-1.5323428) q[2];
sx q[2];
rz(0.65940801) q[2];
rz(-0.13255969) q[3];
sx q[3];
rz(-2.5949251) q[3];
sx q[3];
rz(2.4373655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34683126) q[0];
sx q[0];
rz(-0.95585388) q[0];
sx q[0];
rz(-0.26926789) q[0];
rz(1.5003834) q[1];
sx q[1];
rz(-1.8625448) q[1];
sx q[1];
rz(-2.0620811) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0844042) q[0];
sx q[0];
rz(-1.7499754) q[0];
sx q[0];
rz(-0.65177952) q[0];
rz(2.232021) q[2];
sx q[2];
rz(-2.231763) q[2];
sx q[2];
rz(-1.8125507) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0690382) q[1];
sx q[1];
rz(-2.5423706) q[1];
sx q[1];
rz(-0.95675442) q[1];
rz(0.31793873) q[3];
sx q[3];
rz(-1.7722181) q[3];
sx q[3];
rz(-1.458223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.72011224) q[2];
sx q[2];
rz(-2.5613027) q[2];
sx q[2];
rz(-0.87023467) q[2];
rz(-1.8057711) q[3];
sx q[3];
rz(-1.3355052) q[3];
sx q[3];
rz(-2.465872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0077591) q[0];
sx q[0];
rz(-0.02820153) q[0];
sx q[0];
rz(-0.61122417) q[0];
rz(-2.4080343) q[1];
sx q[1];
rz(-1.4708142) q[1];
sx q[1];
rz(0.38331097) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54827842) q[0];
sx q[0];
rz(-2.3735974) q[0];
sx q[0];
rz(0.40582121) q[0];
x q[1];
rz(-2.3141642) q[2];
sx q[2];
rz(-1.4216004) q[2];
sx q[2];
rz(-2.0713446) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1997879) q[1];
sx q[1];
rz(-2.2082303) q[1];
sx q[1];
rz(0.15934039) q[1];
x q[2];
rz(-2.681045) q[3];
sx q[3];
rz(-0.78118284) q[3];
sx q[3];
rz(-0.51370843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7395301) q[2];
sx q[2];
rz(-0.5548839) q[2];
sx q[2];
rz(-1.0129207) q[2];
rz(0.5976451) q[3];
sx q[3];
rz(-1.4089855) q[3];
sx q[3];
rz(-2.2156318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-1.0585854) q[0];
sx q[0];
rz(-2.7124131) q[0];
sx q[0];
rz(0.035932628) q[0];
rz(1.2220471) q[1];
sx q[1];
rz(-2.4655894) q[1];
sx q[1];
rz(-2.8935208) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1897514) q[0];
sx q[0];
rz(-1.3782256) q[0];
sx q[0];
rz(-1.309881) q[0];
rz(-1.4291968) q[2];
sx q[2];
rz(-1.0842825) q[2];
sx q[2];
rz(1.6213191) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2100468) q[1];
sx q[1];
rz(-1.4074874) q[1];
sx q[1];
rz(1.1985267) q[1];
x q[2];
rz(-1.5733499) q[3];
sx q[3];
rz(-1.4658584) q[3];
sx q[3];
rz(-1.5401287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0849358) q[2];
sx q[2];
rz(-1.0597022) q[2];
sx q[2];
rz(0.53209957) q[2];
rz(0.45446864) q[3];
sx q[3];
rz(-1.5458509) q[3];
sx q[3];
rz(1.2940297) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4249975) q[0];
sx q[0];
rz(-2.0918562) q[0];
sx q[0];
rz(2.9014034) q[0];
rz(2.5364618) q[1];
sx q[1];
rz(-1.3477707) q[1];
sx q[1];
rz(0.64632195) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8511348) q[0];
sx q[0];
rz(-0.86800985) q[0];
sx q[0];
rz(-0.39726195) q[0];
x q[1];
rz(1.4736299) q[2];
sx q[2];
rz(-0.36280131) q[2];
sx q[2];
rz(-1.7513315) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.60257116) q[1];
sx q[1];
rz(-2.638431) q[1];
sx q[1];
rz(2.7064067) q[1];
rz(-pi) q[2];
rz(-1.8490995) q[3];
sx q[3];
rz(-0.67453803) q[3];
sx q[3];
rz(1.9703678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0295082) q[2];
sx q[2];
rz(-1.4936451) q[2];
sx q[2];
rz(3.0599111) q[2];
rz(2.2937842) q[3];
sx q[3];
rz(-0.40926465) q[3];
sx q[3];
rz(-2.9186287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34143701) q[0];
sx q[0];
rz(-2.2969022) q[0];
sx q[0];
rz(-0.89168125) q[0];
rz(1.910123) q[1];
sx q[1];
rz(-0.48858085) q[1];
sx q[1];
rz(-2.5317392) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22287908) q[0];
sx q[0];
rz(-1.0787258) q[0];
sx q[0];
rz(-0.52910536) q[0];
x q[1];
rz(-1.5853204) q[2];
sx q[2];
rz(-2.3554554) q[2];
sx q[2];
rz(-2.5753389) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.6285204) q[1];
sx q[1];
rz(-1.8634463) q[1];
sx q[1];
rz(0.77192941) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2256919) q[3];
sx q[3];
rz(-1.5865579) q[3];
sx q[3];
rz(2.070316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8651198) q[2];
sx q[2];
rz(-2.1742994) q[2];
sx q[2];
rz(2.8709732) q[2];
rz(-2.596415) q[3];
sx q[3];
rz(-1.3278278) q[3];
sx q[3];
rz(-2.8934208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(0.64428627) q[0];
sx q[0];
rz(-1.8338642) q[0];
sx q[0];
rz(2.9551031) q[0];
rz(1.2093774) q[1];
sx q[1];
rz(-1.7861563) q[1];
sx q[1];
rz(-0.019651042) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3075382) q[0];
sx q[0];
rz(-0.77132934) q[0];
sx q[0];
rz(2.176009) q[0];
x q[1];
rz(-2.8803359) q[2];
sx q[2];
rz(-1.3105596) q[2];
sx q[2];
rz(0.90554905) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3001067) q[1];
sx q[1];
rz(-1.7238614) q[1];
sx q[1];
rz(2.7859405) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3778565) q[3];
sx q[3];
rz(-2.0621571) q[3];
sx q[3];
rz(-1.9869252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7095715) q[2];
sx q[2];
rz(-2.186128) q[2];
sx q[2];
rz(-1.552399) q[2];
rz(-3.0324557) q[3];
sx q[3];
rz(-1.8692632) q[3];
sx q[3];
rz(-0.2763589) q[3];
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
rz(-0.5149665) q[0];
sx q[0];
rz(-1.4780541) q[0];
sx q[0];
rz(1.9274101) q[0];
rz(1.7264438) q[1];
sx q[1];
rz(-1.0217923) q[1];
sx q[1];
rz(1.6820977) q[1];
rz(-2.5900415) q[2];
sx q[2];
rz(-2.6307905) q[2];
sx q[2];
rz(-0.62340977) q[2];
rz(0.23701238) q[3];
sx q[3];
rz(-2.7494299) q[3];
sx q[3];
rz(1.3361479) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
