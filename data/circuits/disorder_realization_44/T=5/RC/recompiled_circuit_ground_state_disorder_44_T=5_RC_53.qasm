OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3888336) q[0];
sx q[0];
rz(-1.692481) q[0];
sx q[0];
rz(-1.4504855) q[0];
rz(-5.3095498) q[1];
sx q[1];
rz(7.7205478) q[1];
sx q[1];
rz(7.2024495) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0231173) q[0];
sx q[0];
rz(-0.035041172) q[0];
sx q[0];
rz(-2.6316597) q[0];
x q[1];
rz(-0.64627846) q[2];
sx q[2];
rz(-0.15028223) q[2];
sx q[2];
rz(-0.26963797) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3340993) q[1];
sx q[1];
rz(-1.6624644) q[1];
sx q[1];
rz(2.2794072) q[1];
x q[2];
rz(-1.1000865) q[3];
sx q[3];
rz(-1.19095) q[3];
sx q[3];
rz(2.0947411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3322525) q[2];
sx q[2];
rz(-1.4490178) q[2];
sx q[2];
rz(1.7285041) q[2];
rz(2.938802) q[3];
sx q[3];
rz(-1.3821802) q[3];
sx q[3];
rz(-3.0387759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24421144) q[0];
sx q[0];
rz(-1.0845217) q[0];
sx q[0];
rz(0.71075034) q[0];
rz(2.3731025) q[1];
sx q[1];
rz(-1.0667421) q[1];
sx q[1];
rz(-1.01952) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43895129) q[0];
sx q[0];
rz(-2.0176689) q[0];
sx q[0];
rz(1.5269482) q[0];
x q[1];
rz(-1.4065845) q[2];
sx q[2];
rz(-2.7052042) q[2];
sx q[2];
rz(-1.1498888) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.63377127) q[1];
sx q[1];
rz(-1.941676) q[1];
sx q[1];
rz(-0.033286496) q[1];
rz(1.016576) q[3];
sx q[3];
rz(-1.4999522) q[3];
sx q[3];
rz(1.4431825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3780313) q[2];
sx q[2];
rz(-1.2445933) q[2];
sx q[2];
rz(1.2223318) q[2];
rz(-1.9289121) q[3];
sx q[3];
rz(-1.0168394) q[3];
sx q[3];
rz(-0.71162629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0843622) q[0];
sx q[0];
rz(-0.98452345) q[0];
sx q[0];
rz(1.9236176) q[0];
rz(-0.51586622) q[1];
sx q[1];
rz(-2.5936544) q[1];
sx q[1];
rz(2.2191494) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2005916) q[0];
sx q[0];
rz(-3.0709478) q[0];
sx q[0];
rz(-0.61265041) q[0];
rz(2.5848021) q[2];
sx q[2];
rz(-2.3344731) q[2];
sx q[2];
rz(-3.0592164) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.53239142) q[1];
sx q[1];
rz(-2.2431886) q[1];
sx q[1];
rz(-0.81873853) q[1];
x q[2];
rz(1.7397142) q[3];
sx q[3];
rz(-1.8955821) q[3];
sx q[3];
rz(1.5320154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1232274) q[2];
sx q[2];
rz(-1.0845228) q[2];
sx q[2];
rz(-1.8017192) q[2];
rz(-0.54723048) q[3];
sx q[3];
rz(-2.7180143) q[3];
sx q[3];
rz(-0.71119285) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0860586) q[0];
sx q[0];
rz(-2.9965897) q[0];
sx q[0];
rz(-2.9823629) q[0];
rz(0.010146443) q[1];
sx q[1];
rz(-2.1187014) q[1];
sx q[1];
rz(-2.8841282) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2672294) q[0];
sx q[0];
rz(-0.82634514) q[0];
sx q[0];
rz(0.49212014) q[0];
rz(-pi) q[1];
rz(1.6370378) q[2];
sx q[2];
rz(-0.77755723) q[2];
sx q[2];
rz(-2.24868) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5415693) q[1];
sx q[1];
rz(-2.5065055) q[1];
sx q[1];
rz(-0.16819185) q[1];
rz(2.2063755) q[3];
sx q[3];
rz(-0.87163371) q[3];
sx q[3];
rz(-2.2548667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.7633535) q[2];
sx q[2];
rz(-1.8469609) q[2];
sx q[2];
rz(2.6687458) q[2];
rz(-2.4371448) q[3];
sx q[3];
rz(-1.364578) q[3];
sx q[3];
rz(-0.64594597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5064297) q[0];
sx q[0];
rz(-1.1744873) q[0];
sx q[0];
rz(-3.0761062) q[0];
rz(-0.40924117) q[1];
sx q[1];
rz(-2.0114653) q[1];
sx q[1];
rz(1.7154891) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5201841) q[0];
sx q[0];
rz(-1.6637633) q[0];
sx q[0];
rz(2.9220102) q[0];
rz(-pi) q[1];
x q[1];
rz(1.481856) q[2];
sx q[2];
rz(-1.1834025) q[2];
sx q[2];
rz(-0.26814869) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.1864422) q[1];
sx q[1];
rz(-1.3480061) q[1];
sx q[1];
rz(-2.8563315) q[1];
rz(-pi) q[2];
rz(-2.0992005) q[3];
sx q[3];
rz(-2.8082153) q[3];
sx q[3];
rz(2.7957145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9466729) q[2];
sx q[2];
rz(-1.0871004) q[2];
sx q[2];
rz(-1.8348414) q[2];
rz(1.170018) q[3];
sx q[3];
rz(-0.40478671) q[3];
sx q[3];
rz(0.10725966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4829247) q[0];
sx q[0];
rz(-1.1558477) q[0];
sx q[0];
rz(-0.40147993) q[0];
rz(0.67277706) q[1];
sx q[1];
rz(-2.3428226) q[1];
sx q[1];
rz(2.2875517) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6333165) q[0];
sx q[0];
rz(-2.9190953) q[0];
sx q[0];
rz(-1.359324) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9561192) q[2];
sx q[2];
rz(-1.7181686) q[2];
sx q[2];
rz(2.3285248) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6352977) q[1];
sx q[1];
rz(-2.1599401) q[1];
sx q[1];
rz(1.4511257) q[1];
x q[2];
rz(0.60131945) q[3];
sx q[3];
rz(-2.4375705) q[3];
sx q[3];
rz(-0.46425925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.65471571) q[2];
sx q[2];
rz(-2.2701023) q[2];
sx q[2];
rz(-1.1775449) q[2];
rz(-0.55142895) q[3];
sx q[3];
rz(-1.5588372) q[3];
sx q[3];
rz(0.83561713) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7198782) q[0];
sx q[0];
rz(-0.28110176) q[0];
sx q[0];
rz(1.9792492) q[0];
rz(1.1516736) q[1];
sx q[1];
rz(-1.78777) q[1];
sx q[1];
rz(2.2231359) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8480523) q[0];
sx q[0];
rz(-1.2402724) q[0];
sx q[0];
rz(-2.8577515) q[0];
x q[1];
rz(-1.2612016) q[2];
sx q[2];
rz(-2.4221276) q[2];
sx q[2];
rz(1.3165064) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9388401) q[1];
sx q[1];
rz(-1.563464) q[1];
sx q[1];
rz(-1.5640902) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9711668) q[3];
sx q[3];
rz(-1.9099351) q[3];
sx q[3];
rz(-2.2234774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.034417001) q[2];
sx q[2];
rz(-1.3849266) q[2];
sx q[2];
rz(0.50298634) q[2];
rz(2.3668187) q[3];
sx q[3];
rz(-0.46190244) q[3];
sx q[3];
rz(1.7983961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69304943) q[0];
sx q[0];
rz(-2.9149084) q[0];
sx q[0];
rz(-1.6402798) q[0];
rz(2.8003108) q[1];
sx q[1];
rz(-1.4309859) q[1];
sx q[1];
rz(-1.1192082) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6014746) q[0];
sx q[0];
rz(-2.1261423) q[0];
sx q[0];
rz(1.1309654) q[0];
x q[1];
rz(1.0985435) q[2];
sx q[2];
rz(-0.53016463) q[2];
sx q[2];
rz(-1.7165754) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4367366) q[1];
sx q[1];
rz(-2.4697127) q[1];
sx q[1];
rz(0.19499548) q[1];
rz(-pi) q[2];
rz(-2.8265727) q[3];
sx q[3];
rz(-0.4128939) q[3];
sx q[3];
rz(1.9747495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.98562733) q[2];
sx q[2];
rz(-2.1817744) q[2];
sx q[2];
rz(-2.7723374) q[2];
rz(-2.8333832) q[3];
sx q[3];
rz(-2.1741368) q[3];
sx q[3];
rz(1.6109899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8592598) q[0];
sx q[0];
rz(-1.8349324) q[0];
sx q[0];
rz(-2.7985213) q[0];
rz(1.9725017) q[1];
sx q[1];
rz(-1.7770551) q[1];
sx q[1];
rz(1.7128568) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0051826) q[0];
sx q[0];
rz(-2.6652968) q[0];
sx q[0];
rz(-2.5695557) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3071204) q[2];
sx q[2];
rz(-1.021046) q[2];
sx q[2];
rz(-2.2817734) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1395806) q[1];
sx q[1];
rz(-1.6320458) q[1];
sx q[1];
rz(-1.1948001) q[1];
rz(-pi) q[2];
rz(-2.2637706) q[3];
sx q[3];
rz(-1.1980499) q[3];
sx q[3];
rz(-0.28723785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.94065654) q[2];
sx q[2];
rz(-2.9328465) q[2];
sx q[2];
rz(1.3915871) q[2];
rz(2.7348147) q[3];
sx q[3];
rz(-1.6986366) q[3];
sx q[3];
rz(-0.87882915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10130356) q[0];
sx q[0];
rz(-2.5387796) q[0];
sx q[0];
rz(-1.783675) q[0];
rz(1.2376002) q[1];
sx q[1];
rz(-1.0126746) q[1];
sx q[1];
rz(2.3416669) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7163776) q[0];
sx q[0];
rz(-2.3122462) q[0];
sx q[0];
rz(2.8275376) q[0];
rz(-2.5913057) q[2];
sx q[2];
rz(-1.3105416) q[2];
sx q[2];
rz(-1.257892) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.25617304) q[1];
sx q[1];
rz(-2.2347663) q[1];
sx q[1];
rz(-1.7682942) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1198197) q[3];
sx q[3];
rz(-2.6116237) q[3];
sx q[3];
rz(-1.3254904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.37821975) q[2];
sx q[2];
rz(-1.7174481) q[2];
sx q[2];
rz(-3.0685032) q[2];
rz(2.8857005) q[3];
sx q[3];
rz(-0.81232324) q[3];
sx q[3];
rz(-1.488744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8681317) q[0];
sx q[0];
rz(-1.4556226) q[0];
sx q[0];
rz(1.8709394) q[0];
rz(1.3399667) q[1];
sx q[1];
rz(-1.3506964) q[1];
sx q[1];
rz(0.66257308) q[1];
rz(-1.5688098) q[2];
sx q[2];
rz(-0.86138267) q[2];
sx q[2];
rz(-1.3087261) q[2];
rz(-2.5005199) q[3];
sx q[3];
rz(-2.0430293) q[3];
sx q[3];
rz(0.52386491) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
