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
rz(3.102432) q[0];
sx q[0];
rz(10.089212) q[0];
rz(-0.18222624) q[1];
sx q[1];
rz(4.8290904) q[1];
sx q[1];
rz(11.157142) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8614757) q[0];
sx q[0];
rz(-1.404477) q[0];
sx q[0];
rz(-0.54027277) q[0];
x q[1];
rz(2.6970942) q[2];
sx q[2];
rz(-1.2685565) q[2];
sx q[2];
rz(-1.5950749) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3388901) q[1];
sx q[1];
rz(-0.65734184) q[1];
sx q[1];
rz(-2.606462) q[1];
rz(-pi) q[2];
rz(1.208918) q[3];
sx q[3];
rz(-1.3986949) q[3];
sx q[3];
rz(0.73735305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0613784) q[2];
sx q[2];
rz(-0.55157101) q[2];
sx q[2];
rz(-0.75458327) q[2];
rz(1.6825698) q[3];
sx q[3];
rz(-2.2857234) q[3];
sx q[3];
rz(-1.4065019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.9262806) q[0];
sx q[0];
rz(-0.54406852) q[0];
sx q[0];
rz(-1.1625483) q[0];
rz(-2.4303719) q[1];
sx q[1];
rz(-0.7178719) q[1];
sx q[1];
rz(-0.1444764) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4220336) q[0];
sx q[0];
rz(-1.8633909) q[0];
sx q[0];
rz(1.6574615) q[0];
x q[1];
rz(2.2644193) q[2];
sx q[2];
rz(-2.2899004) q[2];
sx q[2];
rz(0.70110496) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4481758) q[1];
sx q[1];
rz(-0.87478335) q[1];
sx q[1];
rz(-1.0820746) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8322093) q[3];
sx q[3];
rz(-0.9321292) q[3];
sx q[3];
rz(-0.84413278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7084536) q[2];
sx q[2];
rz(-1.3044367) q[2];
sx q[2];
rz(-2.8348095) q[2];
rz(2.3661546) q[3];
sx q[3];
rz(-0.38475761) q[3];
sx q[3];
rz(-0.14498372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6590092) q[0];
sx q[0];
rz(-0.44335303) q[0];
sx q[0];
rz(2.3495667) q[0];
rz(-1.5576942) q[1];
sx q[1];
rz(-1.7894824) q[1];
sx q[1];
rz(-0.99383324) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3344468) q[0];
sx q[0];
rz(-0.69753555) q[0];
sx q[0];
rz(-0.56484449) q[0];
rz(2.5517919) q[2];
sx q[2];
rz(-1.6720534) q[2];
sx q[2];
rz(2.4112005) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6043338) q[1];
sx q[1];
rz(-2.6151665) q[1];
sx q[1];
rz(-0.98186214) q[1];
x q[2];
rz(-1.9478071) q[3];
sx q[3];
rz(-2.6600231) q[3];
sx q[3];
rz(-1.6653614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1899679) q[2];
sx q[2];
rz(-0.8364532) q[2];
sx q[2];
rz(0.17568976) q[2];
rz(0.22819337) q[3];
sx q[3];
rz(-0.56932813) q[3];
sx q[3];
rz(-2.6364251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93155414) q[0];
sx q[0];
rz(-2.2966972) q[0];
sx q[0];
rz(1.5799874) q[0];
rz(-2.2700894) q[1];
sx q[1];
rz(-1.5305287) q[1];
sx q[1];
rz(2.9759488) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7270532) q[0];
sx q[0];
rz(-3.139608) q[0];
sx q[0];
rz(0.63568799) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5890938) q[2];
sx q[2];
rz(-0.4457655) q[2];
sx q[2];
rz(-2.183977) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6871275) q[1];
sx q[1];
rz(-1.4923959) q[1];
sx q[1];
rz(2.3509174) q[1];
rz(-pi) q[2];
rz(2.5609162) q[3];
sx q[3];
rz(-1.8044953) q[3];
sx q[3];
rz(3.0729938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6988301) q[2];
sx q[2];
rz(-2.5593968) q[2];
sx q[2];
rz(2.3637135) q[2];
rz(-2.2569518) q[3];
sx q[3];
rz(-1.9886465) q[3];
sx q[3];
rz(0.9790498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.00064656249) q[0];
sx q[0];
rz(-1.0362754) q[0];
sx q[0];
rz(2.272814) q[0];
rz(-2.2390305) q[1];
sx q[1];
rz(-1.9156009) q[1];
sx q[1];
rz(2.5696519) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40232274) q[0];
sx q[0];
rz(-0.86193854) q[0];
sx q[0];
rz(-0.92415442) q[0];
x q[1];
rz(-2.4221572) q[2];
sx q[2];
rz(-0.63221064) q[2];
sx q[2];
rz(0.84144634) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2635138) q[1];
sx q[1];
rz(-1.1357508) q[1];
sx q[1];
rz(-0.33532354) q[1];
rz(-1.192351) q[3];
sx q[3];
rz(-1.8009543) q[3];
sx q[3];
rz(2.5251021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.14704412) q[2];
sx q[2];
rz(-0.93783718) q[2];
sx q[2];
rz(0.63329548) q[2];
rz(0.58081943) q[3];
sx q[3];
rz(-1.7827026) q[3];
sx q[3];
rz(-2.4766428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13034114) q[0];
sx q[0];
rz(-2.7524152) q[0];
sx q[0];
rz(1.8036386) q[0];
rz(1.411571) q[1];
sx q[1];
rz(-2.7346225) q[1];
sx q[1];
rz(0.79708797) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10182285) q[0];
sx q[0];
rz(-1.8443786) q[0];
sx q[0];
rz(0.16353971) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6930137) q[2];
sx q[2];
rz(-0.56899348) q[2];
sx q[2];
rz(0.48556604) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.89576777) q[1];
sx q[1];
rz(-1.3029769) q[1];
sx q[1];
rz(0.16606776) q[1];
rz(-pi) q[2];
rz(-2.9607746) q[3];
sx q[3];
rz(-1.0581087) q[3];
sx q[3];
rz(0.11845438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7756614) q[2];
sx q[2];
rz(-0.66880995) q[2];
sx q[2];
rz(-2.456341) q[2];
rz(2.8819486) q[3];
sx q[3];
rz(-2.1681867) q[3];
sx q[3];
rz(1.9297622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24945666) q[0];
sx q[0];
rz(-0.15661713) q[0];
sx q[0];
rz(-0.53949612) q[0];
rz(0.19142137) q[1];
sx q[1];
rz(-1.6554662) q[1];
sx q[1];
rz(-2.6991381) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.060613077) q[0];
sx q[0];
rz(-0.0091305841) q[0];
sx q[0];
rz(1.6083482) q[0];
x q[1];
rz(-0.39789756) q[2];
sx q[2];
rz(-2.1566763) q[2];
sx q[2];
rz(0.90998024) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.328622) q[1];
sx q[1];
rz(-0.96862205) q[1];
sx q[1];
rz(0.41090907) q[1];
x q[2];
rz(2.4367853) q[3];
sx q[3];
rz(-2.2205345) q[3];
sx q[3];
rz(0.55834246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.847984) q[2];
sx q[2];
rz(-2.831735) q[2];
sx q[2];
rz(-0.7027182) q[2];
rz(0.27788776) q[3];
sx q[3];
rz(-2.0746456) q[3];
sx q[3];
rz(1.8028629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8948995) q[0];
sx q[0];
rz(-2.2791857) q[0];
sx q[0];
rz(2.605751) q[0];
rz(0.72227532) q[1];
sx q[1];
rz(-1.7012137) q[1];
sx q[1];
rz(2.3586418) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0412126) q[0];
sx q[0];
rz(-1.6579275) q[0];
sx q[0];
rz(-1.210733) q[0];
x q[1];
rz(-0.39733072) q[2];
sx q[2];
rz(-1.518992) q[2];
sx q[2];
rz(-1.0578294) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4346501) q[1];
sx q[1];
rz(-2.1044113) q[1];
sx q[1];
rz(-2.5360819) q[1];
x q[2];
rz(2.1184741) q[3];
sx q[3];
rz(-1.191282) q[3];
sx q[3];
rz(-1.7051769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4202412) q[2];
sx q[2];
rz(-0.46478096) q[2];
sx q[2];
rz(-0.50344023) q[2];
rz(-2.8069046) q[3];
sx q[3];
rz(-1.8036333) q[3];
sx q[3];
rz(0.23505841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47401416) q[0];
sx q[0];
rz(-2.773371) q[0];
sx q[0];
rz(1.0598805) q[0];
rz(2.8267951) q[1];
sx q[1];
rz(-1.7960725) q[1];
sx q[1];
rz(-1.5816636) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2602291) q[0];
sx q[0];
rz(-2.2405528) q[0];
sx q[0];
rz(1.2854281) q[0];
rz(-pi) q[1];
rz(1.4467116) q[2];
sx q[2];
rz(-0.26657924) q[2];
sx q[2];
rz(-2.3125092) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8851354) q[1];
sx q[1];
rz(-0.68084913) q[1];
sx q[1];
rz(-1.225994) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.022059343) q[3];
sx q[3];
rz(-1.1556871) q[3];
sx q[3];
rz(-0.65655113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9554837) q[2];
sx q[2];
rz(-1.9766221) q[2];
sx q[2];
rz(1.0355518) q[2];
rz(-1.1459076) q[3];
sx q[3];
rz(-1.112554) q[3];
sx q[3];
rz(-0.15315332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2027407) q[0];
sx q[0];
rz(-0.11883141) q[0];
sx q[0];
rz(-0.76549292) q[0];
rz(-1.0368404) q[1];
sx q[1];
rz(-1.2988657) q[1];
sx q[1];
rz(0.11229215) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13829198) q[0];
sx q[0];
rz(-1.1213741) q[0];
sx q[0];
rz(-1.9339191) q[0];
x q[1];
rz(2.3621777) q[2];
sx q[2];
rz(-2.0287435) q[2];
sx q[2];
rz(-3.0061257) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4203305) q[1];
sx q[1];
rz(-0.76618505) q[1];
sx q[1];
rz(1.8618733) q[1];
rz(0.046979172) q[3];
sx q[3];
rz(-2.3571797) q[3];
sx q[3];
rz(1.7681153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.74467337) q[2];
sx q[2];
rz(-1.8645218) q[2];
sx q[2];
rz(-0.12167682) q[2];
rz(-2.7820898) q[3];
sx q[3];
rz(-2.4183141) q[3];
sx q[3];
rz(0.017688964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.551238) q[0];
sx q[0];
rz(-1.4689162) q[0];
sx q[0];
rz(1.3015251) q[0];
rz(1.7474668) q[1];
sx q[1];
rz(-1.9448517) q[1];
sx q[1];
rz(1.3762884) q[1];
rz(2.0040705) q[2];
sx q[2];
rz(-2.1480297) q[2];
sx q[2];
rz(0.99110023) q[2];
rz(1.1351552) q[3];
sx q[3];
rz(-1.6899077) q[3];
sx q[3];
rz(-1.6514889) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
