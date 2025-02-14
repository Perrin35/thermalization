OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0661434) q[0];
sx q[0];
rz(-2.0976522) q[0];
sx q[0];
rz(-0.010308417) q[0];
rz(-0.97996867) q[1];
sx q[1];
rz(4.5865321) q[1];
sx q[1];
rz(10.001339) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86069292) q[0];
sx q[0];
rz(-1.4497978) q[0];
sx q[0];
rz(2.8314986) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2498807) q[2];
sx q[2];
rz(-1.0701961) q[2];
sx q[2];
rz(-2.5426189) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.7652055) q[1];
sx q[1];
rz(-0.55843267) q[1];
sx q[1];
rz(0.33513481) q[1];
x q[2];
rz(2.303894) q[3];
sx q[3];
rz(-1.695096) q[3];
sx q[3];
rz(-1.4364786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6171241) q[2];
sx q[2];
rz(-1.4602129) q[2];
sx q[2];
rz(-1.8560393) q[2];
rz(-1.5517976) q[3];
sx q[3];
rz(-1.0714622) q[3];
sx q[3];
rz(-1.0514528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0153506) q[0];
sx q[0];
rz(-1.9168357) q[0];
sx q[0];
rz(-2.8958564) q[0];
rz(-1.0579146) q[1];
sx q[1];
rz(-1.825288) q[1];
sx q[1];
rz(2.8071075) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4870105) q[0];
sx q[0];
rz(-2.062487) q[0];
sx q[0];
rz(1.0500748) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3804383) q[2];
sx q[2];
rz(-1.2710159) q[2];
sx q[2];
rz(1.8331936) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7868488) q[1];
sx q[1];
rz(-1.1965828) q[1];
sx q[1];
rz(2.0207062) q[1];
x q[2];
rz(-1.2859551) q[3];
sx q[3];
rz(-2.1495616) q[3];
sx q[3];
rz(-0.15742403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1841396) q[2];
sx q[2];
rz(-2.2440971) q[2];
sx q[2];
rz(1.755836) q[2];
rz(0.98006788) q[3];
sx q[3];
rz(-0.46531427) q[3];
sx q[3];
rz(0.050962713) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4362713) q[0];
sx q[0];
rz(-0.956981) q[0];
sx q[0];
rz(1.3901688) q[0];
rz(-2.7888489) q[1];
sx q[1];
rz(-0.91068641) q[1];
sx q[1];
rz(-0.62044755) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.016997) q[0];
sx q[0];
rz(-0.11717883) q[0];
sx q[0];
rz(2.4524053) q[0];
rz(-0.98288713) q[2];
sx q[2];
rz(-1.6574142) q[2];
sx q[2];
rz(0.69332214) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7812278) q[1];
sx q[1];
rz(-1.8958127) q[1];
sx q[1];
rz(2.8311391) q[1];
rz(0.87832344) q[3];
sx q[3];
rz(-1.1583503) q[3];
sx q[3];
rz(0.45468047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0188401) q[2];
sx q[2];
rz(-0.41308013) q[2];
sx q[2];
rz(-1.8943465) q[2];
rz(-2.9774169) q[3];
sx q[3];
rz(-1.434606) q[3];
sx q[3];
rz(-2.5119787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4267047) q[0];
sx q[0];
rz(-0.96788228) q[0];
sx q[0];
rz(-1.8545275) q[0];
rz(0.60028589) q[1];
sx q[1];
rz(-1.3515819) q[1];
sx q[1];
rz(2.2379025) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1290734) q[0];
sx q[0];
rz(-0.12618574) q[0];
sx q[0];
rz(-0.55380765) q[0];
rz(-pi) q[1];
rz(1.6389334) q[2];
sx q[2];
rz(-1.1916416) q[2];
sx q[2];
rz(0.013163002) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0541573) q[1];
sx q[1];
rz(-1.8241) q[1];
sx q[1];
rz(-1.3852055) q[1];
rz(-pi) q[2];
rz(0.65151188) q[3];
sx q[3];
rz(-2.6389696) q[3];
sx q[3];
rz(-2.9732957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6220182) q[2];
sx q[2];
rz(-2.3081686) q[2];
sx q[2];
rz(-2.3681417) q[2];
rz(1.8526239) q[3];
sx q[3];
rz(-0.52923146) q[3];
sx q[3];
rz(2.4066063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.89161038) q[0];
sx q[0];
rz(-2.1365428) q[0];
sx q[0];
rz(-0.49215677) q[0];
rz(-2.6026169) q[1];
sx q[1];
rz(-1.4424126) q[1];
sx q[1];
rz(-1.0999365) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5163706) q[0];
sx q[0];
rz(-1.3110135) q[0];
sx q[0];
rz(-2.7894839) q[0];
rz(-pi) q[1];
rz(-0.01077588) q[2];
sx q[2];
rz(-2.1824565) q[2];
sx q[2];
rz(1.3792559) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.392994) q[1];
sx q[1];
rz(-1.6637319) q[1];
sx q[1];
rz(0.8108906) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8995011) q[3];
sx q[3];
rz(-2.0315758) q[3];
sx q[3];
rz(0.24390175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8010572) q[2];
sx q[2];
rz(-0.33581442) q[2];
sx q[2];
rz(2.578793) q[2];
rz(-2.4417012) q[3];
sx q[3];
rz(-1.2908582) q[3];
sx q[3];
rz(0.25115299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7971147) q[0];
sx q[0];
rz(-2.4176702) q[0];
sx q[0];
rz(2.7864454) q[0];
rz(1.9225325) q[1];
sx q[1];
rz(-1.2390169) q[1];
sx q[1];
rz(2.6893137) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66424021) q[0];
sx q[0];
rz(-1.1731804) q[0];
sx q[0];
rz(2.2009322) q[0];
rz(-pi) q[1];
x q[1];
rz(0.38610555) q[2];
sx q[2];
rz(-1.9506644) q[2];
sx q[2];
rz(-2.7136193) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.94231665) q[1];
sx q[1];
rz(-2.8182223) q[1];
sx q[1];
rz(-1.6466696) q[1];
x q[2];
rz(1.4283435) q[3];
sx q[3];
rz(-0.94678426) q[3];
sx q[3];
rz(-2.595866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.63048116) q[2];
sx q[2];
rz(-0.52383542) q[2];
sx q[2];
rz(-3.0044978) q[2];
rz(-1.5591722) q[3];
sx q[3];
rz(-1.1538006) q[3];
sx q[3];
rz(-0.77663511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.1906076) q[0];
sx q[0];
rz(-0.76222104) q[0];
sx q[0];
rz(-1.5451587) q[0];
rz(2.6243788) q[1];
sx q[1];
rz(-1.9904174) q[1];
sx q[1];
rz(0.98701611) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0407437) q[0];
sx q[0];
rz(-0.13253875) q[0];
sx q[0];
rz(0.1610456) q[0];
rz(-1.1140864) q[2];
sx q[2];
rz(-2.4825077) q[2];
sx q[2];
rz(0.50625077) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.92168857) q[1];
sx q[1];
rz(-0.35383407) q[1];
sx q[1];
rz(-2.4262731) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8214885) q[3];
sx q[3];
rz(-2.8134228) q[3];
sx q[3];
rz(1.270297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3551657) q[2];
sx q[2];
rz(-2.1106909) q[2];
sx q[2];
rz(-0.008358566) q[2];
rz(-0.23056325) q[3];
sx q[3];
rz(-1.9356666) q[3];
sx q[3];
rz(-1.6647388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(3.1281328) q[0];
sx q[0];
rz(-3.1210493) q[0];
sx q[0];
rz(1.531456) q[0];
rz(-1.6186591) q[1];
sx q[1];
rz(-1.5459272) q[1];
sx q[1];
rz(1.5628372) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7195936) q[0];
sx q[0];
rz(-0.61924705) q[0];
sx q[0];
rz(-0.25281711) q[0];
rz(-pi) q[1];
rz(-1.7152856) q[2];
sx q[2];
rz(-1.1983216) q[2];
sx q[2];
rz(2.5691751) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.85962379) q[1];
sx q[1];
rz(-2.0323557) q[1];
sx q[1];
rz(2.0412316) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4736797) q[3];
sx q[3];
rz(-1.8258984) q[3];
sx q[3];
rz(0.15984838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.83850399) q[2];
sx q[2];
rz(-1.9946626) q[2];
sx q[2];
rz(0.037503555) q[2];
rz(-2.904902) q[3];
sx q[3];
rz(-0.37875566) q[3];
sx q[3];
rz(-1.8481351) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8732052) q[0];
sx q[0];
rz(-0.78998843) q[0];
sx q[0];
rz(-2.068212) q[0];
rz(2.7997596) q[1];
sx q[1];
rz(-0.57917246) q[1];
sx q[1];
rz(-2.5198708) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9455915) q[0];
sx q[0];
rz(-0.47179963) q[0];
sx q[0];
rz(0.5990754) q[0];
x q[1];
rz(1.4474117) q[2];
sx q[2];
rz(-2.2272791) q[2];
sx q[2];
rz(2.3036164) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.016704233) q[1];
sx q[1];
rz(-2.0604265) q[1];
sx q[1];
rz(-1.917065) q[1];
rz(-pi) q[2];
rz(2.5863566) q[3];
sx q[3];
rz(-2.976417) q[3];
sx q[3];
rz(-3.0495838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6217893) q[2];
sx q[2];
rz(-2.9534464) q[2];
sx q[2];
rz(0.56358799) q[2];
rz(-1.311519) q[3];
sx q[3];
rz(-1.9026285) q[3];
sx q[3];
rz(-0.81609503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1821063) q[0];
sx q[0];
rz(-0.86903787) q[0];
sx q[0];
rz(-2.2667789) q[0];
rz(-1.733571) q[1];
sx q[1];
rz(-0.66649109) q[1];
sx q[1];
rz(-0.63546884) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59171133) q[0];
sx q[0];
rz(-2.3259865) q[0];
sx q[0];
rz(-2.7278103) q[0];
rz(-pi) q[1];
rz(1.9866139) q[2];
sx q[2];
rz(-1.8228616) q[2];
sx q[2];
rz(-0.76463715) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6136352) q[1];
sx q[1];
rz(-1.0571169) q[1];
sx q[1];
rz(-2.6761445) q[1];
x q[2];
rz(-1.0069153) q[3];
sx q[3];
rz(-1.0936519) q[3];
sx q[3];
rz(-1.3729707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1982939) q[2];
sx q[2];
rz(-1.0903) q[2];
sx q[2];
rz(-2.3401006) q[2];
rz(1.21579) q[3];
sx q[3];
rz(-2.2299168) q[3];
sx q[3];
rz(3.099814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53032482) q[0];
sx q[0];
rz(-1.4401191) q[0];
sx q[0];
rz(-2.8339207) q[0];
rz(1.8511741) q[1];
sx q[1];
rz(-2.1967874) q[1];
sx q[1];
rz(-2.8312942) q[1];
rz(-1.4353095) q[2];
sx q[2];
rz(-1.3904962) q[2];
sx q[2];
rz(2.4615859) q[2];
rz(0.17388969) q[3];
sx q[3];
rz(-0.66413838) q[3];
sx q[3];
rz(-0.037353368) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
