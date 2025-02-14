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
rz(2.477159) q[0];
rz(-0.18222624) q[1];
sx q[1];
rz(-1.4540949) q[1];
sx q[1];
rz(1.7323642) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0213703) q[0];
sx q[0];
rz(-0.5628559) q[0];
sx q[0];
rz(-2.8261306) q[0];
rz(-2.514226) q[2];
sx q[2];
rz(-2.6098076) q[2];
sx q[2];
rz(2.5587459) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.69610447) q[1];
sx q[1];
rz(-1.0173807) q[1];
sx q[1];
rz(-1.1958108) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0271985) q[3];
sx q[3];
rz(-0.39908394) q[3];
sx q[3];
rz(1.2582859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0802143) q[2];
sx q[2];
rz(-0.55157101) q[2];
sx q[2];
rz(-0.75458327) q[2];
rz(-1.6825698) q[3];
sx q[3];
rz(-0.85586923) q[3];
sx q[3];
rz(1.7350908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9262806) q[0];
sx q[0];
rz(-2.5975241) q[0];
sx q[0];
rz(1.9790443) q[0];
rz(-2.4303719) q[1];
sx q[1];
rz(-2.4237207) q[1];
sx q[1];
rz(-2.9971163) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1294589) q[0];
sx q[0];
rz(-2.8367865) q[0];
sx q[0];
rz(-2.8617959) q[0];
rz(-0.8501022) q[2];
sx q[2];
rz(-1.0689702) q[2];
sx q[2];
rz(0.36862954) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3857984) q[1];
sx q[1];
rz(-0.82634631) q[1];
sx q[1];
rz(-2.6296294) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30938332) q[3];
sx q[3];
rz(-0.9321292) q[3];
sx q[3];
rz(0.84413278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7084536) q[2];
sx q[2];
rz(-1.3044367) q[2];
sx q[2];
rz(0.30678314) q[2];
rz(-2.3661546) q[3];
sx q[3];
rz(-0.38475761) q[3];
sx q[3];
rz(0.14498372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48258346) q[0];
sx q[0];
rz(-0.44335303) q[0];
sx q[0];
rz(0.79202598) q[0];
rz(-1.5838985) q[1];
sx q[1];
rz(-1.3521103) q[1];
sx q[1];
rz(2.1477594) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4980443) q[0];
sx q[0];
rz(-2.1442765) q[0];
sx q[0];
rz(-1.149096) q[0];
x q[1];
rz(-2.9609072) q[2];
sx q[2];
rz(-0.59741086) q[2];
sx q[2];
rz(-2.1512845) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5842918) q[1];
sx q[1];
rz(-1.2879432) q[1];
sx q[1];
rz(2.0209347) q[1];
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
rz(-1.1899679) q[2];
sx q[2];
rz(-2.3051395) q[2];
sx q[2];
rz(-0.17568976) q[2];
rz(0.22819337) q[3];
sx q[3];
rz(-2.5722645) q[3];
sx q[3];
rz(2.6364251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2100385) q[0];
sx q[0];
rz(-0.84489548) q[0];
sx q[0];
rz(1.5799874) q[0];
rz(-0.87150323) q[1];
sx q[1];
rz(-1.5305287) q[1];
sx q[1];
rz(0.16564381) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77885056) q[0];
sx q[0];
rz(-1.5723933) q[0];
sx q[0];
rz(-1.569618) q[0];
rz(-pi) q[1];
rz(-3.13285) q[2];
sx q[2];
rz(-2.0164818) q[2];
sx q[2];
rz(-2.2042556) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1936744) q[1];
sx q[1];
rz(-0.79371023) q[1];
sx q[1];
rz(3.0315184) q[1];
x q[2];
rz(-1.8481726) q[3];
sx q[3];
rz(-1.0078537) q[3];
sx q[3];
rz(-1.790188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.44276253) q[2];
sx q[2];
rz(-2.5593968) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.00064656249) q[0];
sx q[0];
rz(-2.1053173) q[0];
sx q[0];
rz(2.272814) q[0];
rz(0.90256214) q[1];
sx q[1];
rz(-1.9156009) q[1];
sx q[1];
rz(2.5696519) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8800655) q[0];
sx q[0];
rz(-2.2216317) q[0];
sx q[0];
rz(-2.5291247) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0205021) q[2];
sx q[2];
rz(-1.110198) q[2];
sx q[2];
rz(0.014862343) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2635138) q[1];
sx q[1];
rz(-1.1357508) q[1];
sx q[1];
rz(2.8062691) q[1];
rz(-1.192351) q[3];
sx q[3];
rz(-1.8009543) q[3];
sx q[3];
rz(2.5251021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.14704412) q[2];
sx q[2];
rz(-2.2037555) q[2];
sx q[2];
rz(2.5082972) q[2];
rz(-2.5607732) q[3];
sx q[3];
rz(-1.7827026) q[3];
sx q[3];
rz(0.66494989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0112515) q[0];
sx q[0];
rz(-0.3891775) q[0];
sx q[0];
rz(1.3379541) q[0];
rz(1.411571) q[1];
sx q[1];
rz(-0.4069702) q[1];
sx q[1];
rz(-0.79708797) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4914843) q[0];
sx q[0];
rz(-0.31768018) q[0];
sx q[0];
rz(2.0965212) q[0];
rz(0.077812151) q[2];
sx q[2];
rz(-2.1350265) q[2];
sx q[2];
rz(2.5112453) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.460372) q[1];
sx q[1];
rz(-0.31407323) q[1];
sx q[1];
rz(2.112978) q[1];
x q[2];
rz(-0.18081801) q[3];
sx q[3];
rz(-2.083484) q[3];
sx q[3];
rz(0.11845438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7756614) q[2];
sx q[2];
rz(-2.4727827) q[2];
sx q[2];
rz(2.456341) q[2];
rz(2.8819486) q[3];
sx q[3];
rz(-2.1681867) q[3];
sx q[3];
rz(1.9297622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.892136) q[0];
sx q[0];
rz(-2.9849755) q[0];
sx q[0];
rz(2.6020965) q[0];
rz(2.9501713) q[1];
sx q[1];
rz(-1.4861264) q[1];
sx q[1];
rz(0.44245455) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6689597) q[0];
sx q[0];
rz(-1.5704535) q[0];
sx q[0];
rz(-1.5616722) q[0];
x q[1];
rz(-2.0992958) q[2];
sx q[2];
rz(-2.4467154) q[2];
sx q[2];
rz(-0.26000868) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6574508) q[1];
sx q[1];
rz(-1.2353578) q[1];
sx q[1];
rz(-0.92745933) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4367853) q[3];
sx q[3];
rz(-2.2205345) q[3];
sx q[3];
rz(-0.55834246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.847984) q[2];
sx q[2];
rz(-0.30985761) q[2];
sx q[2];
rz(-2.4388745) q[2];
rz(2.8637049) q[3];
sx q[3];
rz(-2.0746456) q[3];
sx q[3];
rz(1.3387298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8948995) q[0];
sx q[0];
rz(-0.86240697) q[0];
sx q[0];
rz(0.53584164) q[0];
rz(-0.72227532) q[1];
sx q[1];
rz(-1.4403789) q[1];
sx q[1];
rz(2.3586418) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43766681) q[0];
sx q[0];
rz(-1.9294318) q[0];
sx q[0];
rz(-3.0485247) q[0];
x q[1];
rz(-1.5146241) q[2];
sx q[2];
rz(-1.9675641) q[2];
sx q[2];
rz(-0.53469354) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6668876) q[1];
sx q[1];
rz(-2.0829446) q[1];
sx q[1];
rz(-0.94774232) q[1];
x q[2];
rz(0.91714717) q[3];
sx q[3];
rz(-0.65509812) q[3];
sx q[3];
rz(-0.411471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(2.9065342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6675785) q[0];
sx q[0];
rz(-2.773371) q[0];
sx q[0];
rz(-1.0598805) q[0];
rz(2.8267951) q[1];
sx q[1];
rz(-1.3455201) q[1];
sx q[1];
rz(1.5816636) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8813635) q[0];
sx q[0];
rz(-2.2405528) q[0];
sx q[0];
rz(-1.8561646) q[0];
rz(-pi) q[1];
x q[1];
rz(1.306172) q[2];
sx q[2];
rz(-1.5381863) q[2];
sx q[2];
rz(-0.6219686) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0421903) q[1];
sx q[1];
rz(-1.785197) q[1];
sx q[1];
rz(0.91941388) q[1];
rz(-pi) q[2];
rz(-1.6208036) q[3];
sx q[3];
rz(-2.7259318) q[3];
sx q[3];
rz(2.5396944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9554837) q[2];
sx q[2];
rz(-1.9766221) q[2];
sx q[2];
rz(-1.0355518) q[2];
rz(1.995685) q[3];
sx q[3];
rz(-2.0290387) q[3];
sx q[3];
rz(-2.9884393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
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
rz(1.0368404) q[1];
sx q[1];
rz(-1.8427269) q[1];
sx q[1];
rz(0.11229215) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85690391) q[0];
sx q[0];
rz(-0.56988003) q[0];
sx q[0];
rz(2.5068552) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1767581) q[2];
sx q[2];
rz(-2.2529184) q[2];
sx q[2];
rz(1.8471931) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8143554) q[1];
sx q[1];
rz(-0.8443409) q[1];
sx q[1];
rz(-2.8721456) q[1];
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
x q[1];
rz(-0.74467337) q[2];
sx q[2];
rz(-1.2770709) q[2];
sx q[2];
rz(-3.0199158) q[2];
rz(-2.7820898) q[3];
sx q[3];
rz(-0.72327852) q[3];
sx q[3];
rz(3.1239037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59035463) q[0];
sx q[0];
rz(-1.6726765) q[0];
sx q[0];
rz(-1.8400675) q[0];
rz(-1.3941258) q[1];
sx q[1];
rz(-1.9448517) q[1];
sx q[1];
rz(1.3762884) q[1];
rz(-2.0040705) q[2];
sx q[2];
rz(-0.99356298) q[2];
sx q[2];
rz(-2.1504924) q[2];
rz(1.8471424) q[3];
sx q[3];
rz(-0.45062267) q[3];
sx q[3];
rz(-2.9723321) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
