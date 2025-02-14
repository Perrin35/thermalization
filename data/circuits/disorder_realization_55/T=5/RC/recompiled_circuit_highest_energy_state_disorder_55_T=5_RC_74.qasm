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
rz(0.060184181) q[0];
sx q[0];
rz(-2.0826075) q[0];
sx q[0];
rz(-1.0101779) q[0];
rz(-0.8085568) q[1];
sx q[1];
rz(2.8454236) q[1];
sx q[1];
rz(12.256395) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5910019) q[0];
sx q[0];
rz(-0.43568107) q[0];
sx q[0];
rz(0.26623078) q[0];
rz(2.9363382) q[2];
sx q[2];
rz(-1.5927093) q[2];
sx q[2];
rz(0.92682971) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3821564) q[1];
sx q[1];
rz(-1.2431989) q[1];
sx q[1];
rz(2.9345678) q[1];
x q[2];
rz(-2.9397291) q[3];
sx q[3];
rz(-0.84175693) q[3];
sx q[3];
rz(-1.0224467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6628722) q[2];
sx q[2];
rz(-0.24820776) q[2];
sx q[2];
rz(3.0459246) q[2];
rz(3.0293448) q[3];
sx q[3];
rz(-0.87533689) q[3];
sx q[3];
rz(-1.4314502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-1.2410759) q[0];
sx q[0];
rz(-2.2523585) q[0];
sx q[0];
rz(-2.526793) q[0];
rz(-0.88848937) q[1];
sx q[1];
rz(-1.588984) q[1];
sx q[1];
rz(2.6420171) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7730656) q[0];
sx q[0];
rz(-1.8963666) q[0];
sx q[0];
rz(0.21296176) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8728054) q[2];
sx q[2];
rz(-2.6013881) q[2];
sx q[2];
rz(0.42551431) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0608637) q[1];
sx q[1];
rz(-1.7434967) q[1];
sx q[1];
rz(-2.2279146) q[1];
rz(2.8039475) q[3];
sx q[3];
rz(-2.1269375) q[3];
sx q[3];
rz(-0.82586702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8932314) q[2];
sx q[2];
rz(-1.9042559) q[2];
sx q[2];
rz(-0.020817967) q[2];
rz(-1.9013532) q[3];
sx q[3];
rz(-1.6720684) q[3];
sx q[3];
rz(-1.1851236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6388539) q[0];
sx q[0];
rz(-2.8732712) q[0];
sx q[0];
rz(-0.33367208) q[0];
rz(-0.73792136) q[1];
sx q[1];
rz(-1.2260022) q[1];
sx q[1];
rz(-3.0121682) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20724587) q[0];
sx q[0];
rz(-0.6967623) q[0];
sx q[0];
rz(1.8201226) q[0];
rz(-2.4366662) q[2];
sx q[2];
rz(-2.3350596) q[2];
sx q[2];
rz(3.0496917) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8379535) q[1];
sx q[1];
rz(-2.9898414) q[1];
sx q[1];
rz(-2.2084153) q[1];
x q[2];
rz(2.0468726) q[3];
sx q[3];
rz(-0.84730803) q[3];
sx q[3];
rz(-2.0499961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1239329) q[2];
sx q[2];
rz(-1.5256226) q[2];
sx q[2];
rz(1.370149) q[2];
rz(0.74639368) q[3];
sx q[3];
rz(-1.8236225) q[3];
sx q[3];
rz(-3.0588176) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.133404) q[0];
sx q[0];
rz(-1.2097825) q[0];
sx q[0];
rz(-2.5764537) q[0];
rz(-0.42916974) q[1];
sx q[1];
rz(-0.64650911) q[1];
sx q[1];
rz(-0.82829222) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47026004) q[0];
sx q[0];
rz(-1.6506221) q[0];
sx q[0];
rz(-2.3998387) q[0];
rz(-pi) q[1];
rz(-2.8548334) q[2];
sx q[2];
rz(-1.1017403) q[2];
sx q[2];
rz(-3.0931851) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.284621) q[1];
sx q[1];
rz(-1.0880252) q[1];
sx q[1];
rz(1.8716145) q[1];
x q[2];
rz(0.24803646) q[3];
sx q[3];
rz(-1.0674849) q[3];
sx q[3];
rz(2.4865884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4565178) q[2];
sx q[2];
rz(-2.0802616) q[2];
sx q[2];
rz(1.431541) q[2];
rz(2.2020014) q[3];
sx q[3];
rz(-0.51274931) q[3];
sx q[3];
rz(3.140894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13570304) q[0];
sx q[0];
rz(-0.26517427) q[0];
sx q[0];
rz(1.3294719) q[0];
rz(-1.1520518) q[1];
sx q[1];
rz(-1.0877129) q[1];
sx q[1];
rz(-3.1413445) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084293289) q[0];
sx q[0];
rz(-1.3868185) q[0];
sx q[0];
rz(1.3217682) q[0];
rz(0.087281422) q[2];
sx q[2];
rz(-2.1033264) q[2];
sx q[2];
rz(2.0244348) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6387512) q[1];
sx q[1];
rz(-1.0783245) q[1];
sx q[1];
rz(0.3943464) q[1];
x q[2];
rz(1.6432006) q[3];
sx q[3];
rz(-2.3010217) q[3];
sx q[3];
rz(-2.8080468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0033215) q[2];
sx q[2];
rz(-1.891581) q[2];
sx q[2];
rz(3.1357583) q[2];
rz(0.88998574) q[3];
sx q[3];
rz(-0.68166387) q[3];
sx q[3];
rz(-1.1267004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8936845) q[0];
sx q[0];
rz(-1.6981145) q[0];
sx q[0];
rz(-1.4877315) q[0];
rz(-3.1112025) q[1];
sx q[1];
rz(-1.1266212) q[1];
sx q[1];
rz(2.1991275) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1346325) q[0];
sx q[0];
rz(-0.31123268) q[0];
sx q[0];
rz(2.4853021) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29877383) q[2];
sx q[2];
rz(-1.818383) q[2];
sx q[2];
rz(0.58591671) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8928463) q[1];
sx q[1];
rz(-2.5809118) q[1];
sx q[1];
rz(-1.0113202) q[1];
x q[2];
rz(2.0464155) q[3];
sx q[3];
rz(-1.3381357) q[3];
sx q[3];
rz(-1.7847248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5421062) q[2];
sx q[2];
rz(-0.78080559) q[2];
sx q[2];
rz(-1.3055118) q[2];
rz(0.54715884) q[3];
sx q[3];
rz(-1.9360417) q[3];
sx q[3];
rz(0.80593306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.9611573) q[0];
sx q[0];
rz(-1.0338217) q[0];
sx q[0];
rz(-3.0410774) q[0];
rz(0.48249498) q[1];
sx q[1];
rz(-1.1588187) q[1];
sx q[1];
rz(0.85711342) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0470951) q[0];
sx q[0];
rz(-2.408354) q[0];
sx q[0];
rz(-0.49219699) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6202507) q[2];
sx q[2];
rz(-2.2170545) q[2];
sx q[2];
rz(2.6924999) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.65322058) q[1];
sx q[1];
rz(-1.3626507) q[1];
sx q[1];
rz(-0.4076165) q[1];
rz(-pi) q[2];
rz(1.1825652) q[3];
sx q[3];
rz(-0.60743466) q[3];
sx q[3];
rz(-1.5608112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7414005) q[2];
sx q[2];
rz(-2.1705748) q[2];
sx q[2];
rz(2.8391489) q[2];
rz(-0.97964573) q[3];
sx q[3];
rz(-1.0023578) q[3];
sx q[3];
rz(-1.8625331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7804467) q[0];
sx q[0];
rz(-3.0028711) q[0];
sx q[0];
rz(-0.030990344) q[0];
rz(-2.567645) q[1];
sx q[1];
rz(-1.6684883) q[1];
sx q[1];
rz(-1.8642289) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81640667) q[0];
sx q[0];
rz(-1.2379097) q[0];
sx q[0];
rz(-1.9247967) q[0];
rz(-2.8988016) q[2];
sx q[2];
rz(-1.0848248) q[2];
sx q[2];
rz(-2.9832341) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.053169202) q[1];
sx q[1];
rz(-0.49234566) q[1];
sx q[1];
rz(-1.7264992) q[1];
rz(1.7788309) q[3];
sx q[3];
rz(-0.29104189) q[3];
sx q[3];
rz(-2.7179949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.10451) q[2];
sx q[2];
rz(-0.13585486) q[2];
sx q[2];
rz(0.48745298) q[2];
rz(2.4449352) q[3];
sx q[3];
rz(-0.928855) q[3];
sx q[3];
rz(3.0776183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0697407) q[0];
sx q[0];
rz(-0.47436473) q[0];
sx q[0];
rz(-2.790614) q[0];
rz(-0.17414302) q[1];
sx q[1];
rz(-1.6330279) q[1];
sx q[1];
rz(-1.7399656) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.021371776) q[0];
sx q[0];
rz(-1.3231771) q[0];
sx q[0];
rz(-1.6970474) q[0];
rz(-pi) q[1];
rz(2.0171595) q[2];
sx q[2];
rz(-0.87136641) q[2];
sx q[2];
rz(1.0126737) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.40608866) q[1];
sx q[1];
rz(-1.5998518) q[1];
sx q[1];
rz(-1.7543704) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.97934874) q[3];
sx q[3];
rz(-0.92100793) q[3];
sx q[3];
rz(1.5486919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.031124) q[2];
sx q[2];
rz(-0.68507552) q[2];
sx q[2];
rz(2.5049211) q[2];
rz(-2.0311671) q[3];
sx q[3];
rz(-1.5444376) q[3];
sx q[3];
rz(-0.5184263) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3670032) q[0];
sx q[0];
rz(-2.84802) q[0];
sx q[0];
rz(-1.0585744) q[0];
rz(0.038657945) q[1];
sx q[1];
rz(-1.5267173) q[1];
sx q[1];
rz(-1.0640594) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3389726) q[0];
sx q[0];
rz(-1.5741328) q[0];
sx q[0];
rz(0.0044857684) q[0];
rz(-pi) q[1];
rz(-2.8644807) q[2];
sx q[2];
rz(-1.849035) q[2];
sx q[2];
rz(0.81467512) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7239889) q[1];
sx q[1];
rz(-3.0688802) q[1];
sx q[1];
rz(-0.59007187) q[1];
rz(-pi) q[2];
rz(0.070399447) q[3];
sx q[3];
rz(-0.46050763) q[3];
sx q[3];
rz(1.2706336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.71558636) q[2];
sx q[2];
rz(-2.6435489) q[2];
sx q[2];
rz(3.0675724) q[2];
rz(2.7600539) q[3];
sx q[3];
rz(-1.8266725) q[3];
sx q[3];
rz(0.35416245) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1114125) q[0];
sx q[0];
rz(-1.3127865) q[0];
sx q[0];
rz(2.4319613) q[0];
rz(-2.5297655) q[1];
sx q[1];
rz(-2.430293) q[1];
sx q[1];
rz(-1.5878955) q[1];
rz(-0.62803531) q[2];
sx q[2];
rz(-2.5428094) q[2];
sx q[2];
rz(1.6369292) q[2];
rz(-0.5027708) q[3];
sx q[3];
rz(-2.3270861) q[3];
sx q[3];
rz(2.6181639) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
