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
rz(-1.7323642) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7519386) q[0];
sx q[0];
rz(-1.03878) q[0];
sx q[0];
rz(-1.3774894) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44449842) q[2];
sx q[2];
rz(-1.2685565) q[2];
sx q[2];
rz(1.5465178) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4709027) q[1];
sx q[1];
rz(-1.2539314) q[1];
sx q[1];
rz(-0.58611996) q[1];
rz(-2.9578311) q[3];
sx q[3];
rz(-1.9270883) q[3];
sx q[3];
rz(0.7687062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0802143) q[2];
sx q[2];
rz(-0.55157101) q[2];
sx q[2];
rz(0.75458327) q[2];
rz(-1.6825698) q[3];
sx q[3];
rz(-2.2857234) q[3];
sx q[3];
rz(1.4065019) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.215312) q[0];
sx q[0];
rz(-0.54406852) q[0];
sx q[0];
rz(1.9790443) q[0];
rz(-2.4303719) q[1];
sx q[1];
rz(-0.7178719) q[1];
sx q[1];
rz(-0.1444764) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9677749) q[0];
sx q[0];
rz(-1.6537695) q[0];
sx q[0];
rz(2.847958) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2644193) q[2];
sx q[2];
rz(-0.85169221) q[2];
sx q[2];
rz(-2.4404877) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6934168) q[1];
sx q[1];
rz(-2.2668093) q[1];
sx q[1];
rz(1.0820746) q[1];
rz(-pi) q[2];
rz(0.90869168) q[3];
sx q[3];
rz(-1.8177563) q[3];
sx q[3];
rz(-2.2266362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7084536) q[2];
sx q[2];
rz(-1.837156) q[2];
sx q[2];
rz(-2.8348095) q[2];
rz(-0.77543801) q[3];
sx q[3];
rz(-2.756835) q[3];
sx q[3];
rz(-2.9966089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68849148) q[0];
sx q[0];
rz(-1.2198041) q[0];
sx q[0];
rz(2.5255802) q[0];
x q[1];
rz(0.58980073) q[2];
sx q[2];
rz(-1.4695393) q[2];
sx q[2];
rz(2.4112005) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6043338) q[1];
sx q[1];
rz(-0.52642614) q[1];
sx q[1];
rz(2.1597305) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9478071) q[3];
sx q[3];
rz(-0.48156958) q[3];
sx q[3];
rz(1.4762312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9516248) q[2];
sx q[2];
rz(-2.3051395) q[2];
sx q[2];
rz(2.9659029) q[2];
rz(0.22819337) q[3];
sx q[3];
rz(-0.56932813) q[3];
sx q[3];
rz(-2.6364251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2100385) q[0];
sx q[0];
rz(-0.84489548) q[0];
sx q[0];
rz(1.5799874) q[0];
rz(0.87150323) q[1];
sx q[1];
rz(-1.5305287) q[1];
sx q[1];
rz(-0.16564381) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79194389) q[0];
sx q[0];
rz(-1.5719746) q[0];
sx q[0];
rz(0.0015969413) q[0];
rz(0.0087426337) q[2];
sx q[2];
rz(-2.0164818) q[2];
sx q[2];
rz(-2.2042556) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.45446515) q[1];
sx q[1];
rz(-1.6491967) q[1];
sx q[1];
rz(2.3509174) q[1];
rz(-pi) q[2];
rz(-0.58067643) q[3];
sx q[3];
rz(-1.3370974) q[3];
sx q[3];
rz(-3.0729938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6988301) q[2];
sx q[2];
rz(-2.5593968) q[2];
sx q[2];
rz(-2.3637135) q[2];
rz(-0.88464087) q[3];
sx q[3];
rz(-1.1529461) q[3];
sx q[3];
rz(-2.1625429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
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
rz(0.90256214) q[1];
sx q[1];
rz(-1.2259918) q[1];
sx q[1];
rz(-2.5696519) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40232274) q[0];
sx q[0];
rz(-0.86193854) q[0];
sx q[0];
rz(2.2174382) q[0];
rz(-0.71943546) q[2];
sx q[2];
rz(-0.63221064) q[2];
sx q[2];
rz(-0.84144634) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8780788) q[1];
sx q[1];
rz(-1.1357508) q[1];
sx q[1];
rz(0.33532354) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24700332) q[3];
sx q[3];
rz(-1.2028143) q[3];
sx q[3];
rz(2.0968269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9945485) q[2];
sx q[2];
rz(-0.93783718) q[2];
sx q[2];
rz(2.5082972) q[2];
rz(-0.58081943) q[3];
sx q[3];
rz(-1.3588901) q[3];
sx q[3];
rz(-2.4766428) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0112515) q[0];
sx q[0];
rz(-2.7524152) q[0];
sx q[0];
rz(-1.3379541) q[0];
rz(-1.411571) q[1];
sx q[1];
rz(-0.4069702) q[1];
sx q[1];
rz(-2.3445047) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4914843) q[0];
sx q[0];
rz(-0.31768018) q[0];
sx q[0];
rz(-1.0450715) q[0];
x q[1];
rz(0.077812151) q[2];
sx q[2];
rz(-1.0065662) q[2];
sx q[2];
rz(0.63034731) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.460372) q[1];
sx q[1];
rz(-2.8275194) q[1];
sx q[1];
rz(-1.0286147) q[1];
rz(-pi) q[2];
x q[2];
rz(0.18081801) q[3];
sx q[3];
rz(-2.083484) q[3];
sx q[3];
rz(3.0231383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.36593124) q[2];
sx q[2];
rz(-0.66880995) q[2];
sx q[2];
rz(0.68525165) q[2];
rz(0.25964409) q[3];
sx q[3];
rz(-2.1681867) q[3];
sx q[3];
rz(-1.9297622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.098166532) q[0];
sx q[0];
rz(-1.5799205) q[0];
sx q[0];
rz(0.00034279963) q[0];
x q[1];
rz(-2.0992958) q[2];
sx q[2];
rz(-0.69487725) q[2];
sx q[2];
rz(0.26000868) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6574508) q[1];
sx q[1];
rz(-1.9062348) q[1];
sx q[1];
rz(2.2141333) q[1];
rz(2.4367853) q[3];
sx q[3];
rz(-2.2205345) q[3];
sx q[3];
rz(0.55834246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29360867) q[2];
sx q[2];
rz(-2.831735) q[2];
sx q[2];
rz(0.7027182) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8948995) q[0];
sx q[0];
rz(-0.86240697) q[0];
sx q[0];
rz(-2.605751) q[0];
rz(2.4193173) q[1];
sx q[1];
rz(-1.7012137) q[1];
sx q[1];
rz(0.78295082) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0412126) q[0];
sx q[0];
rz(-1.6579275) q[0];
sx q[0];
rz(1.9308596) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0083904) q[2];
sx q[2];
rz(-0.40051546) q[2];
sx q[2];
rz(-0.39019624) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6668876) q[1];
sx q[1];
rz(-2.0829446) q[1];
sx q[1];
rz(0.94774232) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.43704982) q[3];
sx q[3];
rz(-2.0756222) q[3];
sx q[3];
rz(0.35660353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7213514) q[2];
sx q[2];
rz(-2.6768117) q[2];
sx q[2];
rz(-2.6381524) q[2];
rz(2.8069046) q[3];
sx q[3];
rz(-1.8036333) q[3];
sx q[3];
rz(-0.23505841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
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
rz(-1.3455201) q[1];
sx q[1];
rz(1.559929) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2602291) q[0];
sx q[0];
rz(-0.90103982) q[0];
sx q[0];
rz(-1.2854281) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.694881) q[2];
sx q[2];
rz(-0.26657924) q[2];
sx q[2];
rz(0.82908344) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.68941249) q[1];
sx q[1];
rz(-0.93675568) q[1];
sx q[1];
rz(-2.874332) q[1];
rz(-pi) q[2];
rz(1.1555973) q[3];
sx q[3];
rz(-1.590982) q[3];
sx q[3];
rz(0.92314271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9554837) q[2];
sx q[2];
rz(-1.1649705) q[2];
sx q[2];
rz(1.0355518) q[2];
rz(1.995685) q[3];
sx q[3];
rz(-1.112554) q[3];
sx q[3];
rz(-0.15315332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2027407) q[0];
sx q[0];
rz(-3.0227612) q[0];
sx q[0];
rz(0.76549292) q[0];
rz(-1.0368404) q[1];
sx q[1];
rz(-1.2988657) q[1];
sx q[1];
rz(0.11229215) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.87269) q[0];
sx q[0];
rz(-1.8964452) q[0];
sx q[0];
rz(-2.6652357) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61157558) q[2];
sx q[2];
rz(-2.2627352) q[2];
sx q[2];
rz(1.0146146) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.3272373) q[1];
sx q[1];
rz(-0.8443409) q[1];
sx q[1];
rz(-2.8721456) q[1];
rz(-pi) q[2];
rz(1.6176315) q[3];
sx q[3];
rz(-2.3541087) q[3];
sx q[3];
rz(1.3071527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.74467337) q[2];
sx q[2];
rz(-1.2770709) q[2];
sx q[2];
rz(0.12167682) q[2];
rz(-0.35950288) q[3];
sx q[3];
rz(-2.4183141) q[3];
sx q[3];
rz(-0.017688964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.551238) q[0];
sx q[0];
rz(-1.6726765) q[0];
sx q[0];
rz(-1.8400675) q[0];
rz(1.3941258) q[1];
sx q[1];
rz(-1.196741) q[1];
sx q[1];
rz(-1.7653042) q[1];
rz(-0.57264282) q[2];
sx q[2];
rz(-2.4349458) q[2];
sx q[2];
rz(-1.4473421) q[2];
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
