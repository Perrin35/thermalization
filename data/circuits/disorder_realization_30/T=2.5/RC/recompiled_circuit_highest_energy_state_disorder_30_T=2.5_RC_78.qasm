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
rz(-2.8673515) q[0];
sx q[0];
rz(-0.40617901) q[0];
sx q[0];
rz(2.9131373) q[0];
rz(0.38415456) q[1];
sx q[1];
rz(-1.471712) q[1];
sx q[1];
rz(0.44878146) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81918994) q[0];
sx q[0];
rz(-2.7459641) q[0];
sx q[0];
rz(-0.36901335) q[0];
rz(3.0182462) q[2];
sx q[2];
rz(-1.1575067) q[2];
sx q[2];
rz(1.1878428) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8463414) q[1];
sx q[1];
rz(-1.59473) q[1];
sx q[1];
rz(-1.5071391) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5767831) q[3];
sx q[3];
rz(-1.5428169) q[3];
sx q[3];
rz(0.81708032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7479129) q[2];
sx q[2];
rz(-0.23466514) q[2];
sx q[2];
rz(0.62825424) q[2];
rz(1.0660394) q[3];
sx q[3];
rz(-0.99848905) q[3];
sx q[3];
rz(1.1294686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5265441) q[0];
sx q[0];
rz(-0.55416179) q[0];
sx q[0];
rz(2.2224485) q[0];
rz(2.9650086) q[1];
sx q[1];
rz(-1.44839) q[1];
sx q[1];
rz(-0.44763705) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31792163) q[0];
sx q[0];
rz(-1.5632251) q[0];
sx q[0];
rz(-1.5764144) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7451374) q[2];
sx q[2];
rz(-1.8282991) q[2];
sx q[2];
rz(-2.5776209) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0470878) q[1];
sx q[1];
rz(-0.89549235) q[1];
sx q[1];
rz(1.7750791) q[1];
rz(-2.6352386) q[3];
sx q[3];
rz(-1.10276) q[3];
sx q[3];
rz(-2.9705257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.24719079) q[2];
sx q[2];
rz(-1.2195419) q[2];
sx q[2];
rz(2.2578237) q[2];
rz(-0.73924685) q[3];
sx q[3];
rz(-1.2284307) q[3];
sx q[3];
rz(1.242274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3833375) q[0];
sx q[0];
rz(-1.0403591) q[0];
sx q[0];
rz(0.0058767949) q[0];
rz(-0.30644304) q[1];
sx q[1];
rz(-1.076661) q[1];
sx q[1];
rz(1.2642911) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5976083) q[0];
sx q[0];
rz(-0.48564821) q[0];
sx q[0];
rz(0.66868725) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1074171) q[2];
sx q[2];
rz(-1.4219936) q[2];
sx q[2];
rz(0.11916313) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2325461) q[1];
sx q[1];
rz(-1.0139483) q[1];
sx q[1];
rz(1.805748) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39671582) q[3];
sx q[3];
rz(-1.0457888) q[3];
sx q[3];
rz(-2.5535312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3197202) q[2];
sx q[2];
rz(-1.7629273) q[2];
sx q[2];
rz(2.482448) q[2];
rz(1.871073) q[3];
sx q[3];
rz(-2.8489385) q[3];
sx q[3];
rz(1.6856153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(1.4996516) q[0];
sx q[0];
rz(-0.15436509) q[0];
sx q[0];
rz(-0.9374215) q[0];
rz(1.5501267) q[1];
sx q[1];
rz(-0.66842404) q[1];
sx q[1];
rz(-0.74916565) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8892775) q[0];
sx q[0];
rz(-1.2752644) q[0];
sx q[0];
rz(-0.79177971) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90355166) q[2];
sx q[2];
rz(-1.0507116) q[2];
sx q[2];
rz(1.5423216) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8038102) q[1];
sx q[1];
rz(-0.69548464) q[1];
sx q[1];
rz(-1.2616322) q[1];
rz(-0.22386287) q[3];
sx q[3];
rz(-0.77896691) q[3];
sx q[3];
rz(1.9958441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.59430993) q[2];
sx q[2];
rz(-0.82188934) q[2];
sx q[2];
rz(-1.8853356) q[2];
rz(-0.78993434) q[3];
sx q[3];
rz(-2.200685) q[3];
sx q[3];
rz(-1.1676211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.054852) q[0];
sx q[0];
rz(-0.60844839) q[0];
sx q[0];
rz(-0.21859455) q[0];
rz(-0.78544468) q[1];
sx q[1];
rz(-2.0365069) q[1];
sx q[1];
rz(1.2931664) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47951133) q[0];
sx q[0];
rz(-1.1089227) q[0];
sx q[0];
rz(0.9734459) q[0];
x q[1];
rz(-0.082581994) q[2];
sx q[2];
rz(-1.3090927) q[2];
sx q[2];
rz(-2.341389) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.47673273) q[1];
sx q[1];
rz(-0.97366758) q[1];
sx q[1];
rz(2.4684811) q[1];
x q[2];
rz(-0.038083548) q[3];
sx q[3];
rz(-1.3553737) q[3];
sx q[3];
rz(2.0652536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.834632) q[2];
sx q[2];
rz(-2.0768276) q[2];
sx q[2];
rz(1.628423) q[2];
rz(-2.4026134) q[3];
sx q[3];
rz(-1.8578953) q[3];
sx q[3];
rz(1.4449323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4479248) q[0];
sx q[0];
rz(-1.6532927) q[0];
sx q[0];
rz(0.36201763) q[0];
rz(-1.1297049) q[1];
sx q[1];
rz(-1.2672707) q[1];
sx q[1];
rz(2.3806908) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7603455) q[0];
sx q[0];
rz(-1.3804304) q[0];
sx q[0];
rz(-2.3668853) q[0];
rz(-0.82693066) q[2];
sx q[2];
rz(-1.880135) q[2];
sx q[2];
rz(0.29238809) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7089272) q[1];
sx q[1];
rz(-0.58307465) q[1];
sx q[1];
rz(1.9645755) q[1];
x q[2];
rz(1.170916) q[3];
sx q[3];
rz(-0.99046773) q[3];
sx q[3];
rz(2.1023242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0249944) q[2];
sx q[2];
rz(-2.8572539) q[2];
sx q[2];
rz(1.9653758) q[2];
rz(-2.1182012) q[3];
sx q[3];
rz(-2.0659645) q[3];
sx q[3];
rz(0.37337506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0244331) q[0];
sx q[0];
rz(-0.42851448) q[0];
sx q[0];
rz(-1.3065216) q[0];
rz(0.10082968) q[1];
sx q[1];
rz(-1.3918624) q[1];
sx q[1];
rz(0.94591013) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8656509) q[0];
sx q[0];
rz(-1.8748594) q[0];
sx q[0];
rz(-1.434554) q[0];
rz(-pi) q[1];
rz(-0.77581279) q[2];
sx q[2];
rz(-2.472942) q[2];
sx q[2];
rz(2.5324627) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1506699) q[1];
sx q[1];
rz(-1.4353818) q[1];
sx q[1];
rz(-1.6305429) q[1];
x q[2];
rz(-0.85695388) q[3];
sx q[3];
rz(-1.3009239) q[3];
sx q[3];
rz(1.2601579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0239608) q[2];
sx q[2];
rz(-1.5317711) q[2];
sx q[2];
rz(-2.7839933) q[2];
rz(3.0910953) q[3];
sx q[3];
rz(-0.42323083) q[3];
sx q[3];
rz(0.38664204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9015775) q[0];
sx q[0];
rz(-1.3192663) q[0];
sx q[0];
rz(-2.5460119) q[0];
rz(1.2254084) q[1];
sx q[1];
rz(-1.775368) q[1];
sx q[1];
rz(0.43232408) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8416467) q[0];
sx q[0];
rz(-0.52021813) q[0];
sx q[0];
rz(-2.8033517) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.064266845) q[2];
sx q[2];
rz(-1.7241577) q[2];
sx q[2];
rz(-2.4118285) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3485175) q[1];
sx q[1];
rz(-1.5695453) q[1];
sx q[1];
rz(-2.5743471) q[1];
x q[2];
rz(-2.0995977) q[3];
sx q[3];
rz(-2.3862958) q[3];
sx q[3];
rz(-0.55213682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3782184) q[2];
sx q[2];
rz(-1.8209527) q[2];
sx q[2];
rz(2.1136368) q[2];
rz(1.4823191) q[3];
sx q[3];
rz(-1.4488528) q[3];
sx q[3];
rz(-0.75215522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1886275) q[0];
sx q[0];
rz(-1.9099706) q[0];
sx q[0];
rz(-2.6787483) q[0];
rz(-1.9075182) q[1];
sx q[1];
rz(-1.2214829) q[1];
sx q[1];
rz(5*pi/8) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6825487) q[0];
sx q[0];
rz(-0.17691806) q[0];
sx q[0];
rz(1.4260308) q[0];
x q[1];
rz(0.79754957) q[2];
sx q[2];
rz(-1.2204959) q[2];
sx q[2];
rz(2.1458117) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6057889) q[1];
sx q[1];
rz(-2.2767496) q[1];
sx q[1];
rz(-0.16210254) q[1];
rz(-pi) q[2];
x q[2];
rz(0.34171243) q[3];
sx q[3];
rz(-1.4572506) q[3];
sx q[3];
rz(0.36817238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6356293) q[2];
sx q[2];
rz(-1.7545981) q[2];
sx q[2];
rz(-0.76816922) q[2];
rz(-2.4913037) q[3];
sx q[3];
rz(-2.8322329) q[3];
sx q[3];
rz(-0.55997854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3637417) q[0];
sx q[0];
rz(-0.44054458) q[0];
sx q[0];
rz(-1.1458696) q[0];
rz(-0.56529415) q[1];
sx q[1];
rz(-2.4259613) q[1];
sx q[1];
rz(2.735945) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.250428) q[0];
sx q[0];
rz(-1.2235723) q[0];
sx q[0];
rz(-0.77607147) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1937386) q[2];
sx q[2];
rz(-1.8296281) q[2];
sx q[2];
rz(1.7155312) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1310869) q[1];
sx q[1];
rz(-1.0684135) q[1];
sx q[1];
rz(1.2817863) q[1];
x q[2];
rz(-1.6784366) q[3];
sx q[3];
rz(-1.1336796) q[3];
sx q[3];
rz(2.6588065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.046935) q[2];
sx q[2];
rz(-1.2311225) q[2];
sx q[2];
rz(-0.51935736) q[2];
rz(2.6804067) q[3];
sx q[3];
rz(-1.0591732) q[3];
sx q[3];
rz(-0.69380277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2431348) q[0];
sx q[0];
rz(-1.5903789) q[0];
sx q[0];
rz(-0.90573885) q[0];
rz(1.4836736) q[1];
sx q[1];
rz(-1.1653733) q[1];
sx q[1];
rz(1.9001874) q[1];
rz(2.1853191) q[2];
sx q[2];
rz(-1.6884138) q[2];
sx q[2];
rz(-1.5547167) q[2];
rz(0.46403454) q[3];
sx q[3];
rz(-1.1190718) q[3];
sx q[3];
rz(-3.0116871) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
