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
rz(0.54654044) q[0];
sx q[0];
rz(-2.9822783) q[0];
sx q[0];
rz(1.8848609) q[0];
rz(-2.9873084) q[1];
sx q[1];
rz(-1.0839387) q[1];
sx q[1];
rz(1.7854324) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0935554) q[0];
sx q[0];
rz(-1.1021656) q[0];
sx q[0];
rz(-0.2524391) q[0];
x q[1];
rz(-1.5370813) q[2];
sx q[2];
rz(-1.7236934) q[2];
sx q[2];
rz(-0.70399415) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6068597) q[1];
sx q[1];
rz(-2.4280425) q[1];
sx q[1];
rz(-0.092852863) q[1];
rz(-pi) q[2];
rz(0.7103997) q[3];
sx q[3];
rz(-2.5419639) q[3];
sx q[3];
rz(2.8192632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.95947444) q[2];
sx q[2];
rz(-1.139816) q[2];
sx q[2];
rz(-0.095890447) q[2];
rz(0.31283665) q[3];
sx q[3];
rz(-1.0597119) q[3];
sx q[3];
rz(-0.50725168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7392015) q[0];
sx q[0];
rz(-1.6283789) q[0];
sx q[0];
rz(1.8550523) q[0];
rz(-1.3013499) q[1];
sx q[1];
rz(-1.8845314) q[1];
sx q[1];
rz(1.7110862) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.032784) q[0];
sx q[0];
rz(-2.6258402) q[0];
sx q[0];
rz(1.4620251) q[0];
rz(-pi) q[1];
x q[1];
rz(0.44746676) q[2];
sx q[2];
rz(-1.1612084) q[2];
sx q[2];
rz(-1.8597395) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.41005653) q[1];
sx q[1];
rz(-1.6997196) q[1];
sx q[1];
rz(1.6881264) q[1];
rz(-pi) q[2];
rz(-1.3271403) q[3];
sx q[3];
rz(-2.2639963) q[3];
sx q[3];
rz(1.6125295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3133189) q[2];
sx q[2];
rz(-2.2602849) q[2];
sx q[2];
rz(2.4271097) q[2];
rz(-1.03164) q[3];
sx q[3];
rz(-0.98543006) q[3];
sx q[3];
rz(0.45903444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76630509) q[0];
sx q[0];
rz(-2.3630688) q[0];
sx q[0];
rz(0.2226204) q[0];
rz(-0.34046945) q[1];
sx q[1];
rz(-1.4687931) q[1];
sx q[1];
rz(2.8276843) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8720377) q[0];
sx q[0];
rz(-1.1133218) q[0];
sx q[0];
rz(0.78367718) q[0];
rz(-pi) q[1];
x q[1];
rz(0.48431335) q[2];
sx q[2];
rz(-0.90171684) q[2];
sx q[2];
rz(2.1824257) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7463689) q[1];
sx q[1];
rz(-0.88942553) q[1];
sx q[1];
rz(1.1716151) q[1];
rz(-pi) q[2];
rz(-2.7438358) q[3];
sx q[3];
rz(-0.53103775) q[3];
sx q[3];
rz(-0.49285313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8935304) q[2];
sx q[2];
rz(-1.6174822) q[2];
sx q[2];
rz(-1.0244055) q[2];
rz(-1.8223358) q[3];
sx q[3];
rz(-0.28194591) q[3];
sx q[3];
rz(-2.8446741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1612741) q[0];
sx q[0];
rz(-1.4264822) q[0];
sx q[0];
rz(-1.0262161) q[0];
rz(-0.16464344) q[1];
sx q[1];
rz(-2.7963729) q[1];
sx q[1];
rz(2.3688597) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5948295) q[0];
sx q[0];
rz(-1.9332316) q[0];
sx q[0];
rz(-2.4700463) q[0];
rz(0.35773678) q[2];
sx q[2];
rz(-1.5669723) q[2];
sx q[2];
rz(1.2847023) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0028006) q[1];
sx q[1];
rz(-1.7874663) q[1];
sx q[1];
rz(-0.8817706) q[1];
rz(1.5192612) q[3];
sx q[3];
rz(-1.5665352) q[3];
sx q[3];
rz(2.2770027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.905895) q[2];
sx q[2];
rz(-2.5696281) q[2];
sx q[2];
rz(0.31615654) q[2];
rz(0.21137992) q[3];
sx q[3];
rz(-1.6930765) q[3];
sx q[3];
rz(-1.0630509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0036156) q[0];
sx q[0];
rz(-1.6162385) q[0];
sx q[0];
rz(0.77950087) q[0];
rz(1.4315073) q[1];
sx q[1];
rz(-1.5305488) q[1];
sx q[1];
rz(0.12758189) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9923926) q[0];
sx q[0];
rz(-3.0221327) q[0];
sx q[0];
rz(-1.275927) q[0];
x q[1];
rz(-0.76697369) q[2];
sx q[2];
rz(-0.80000814) q[2];
sx q[2];
rz(-2.6376079) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7320805) q[1];
sx q[1];
rz(-1.8276375) q[1];
sx q[1];
rz(-1.9496296) q[1];
rz(-pi) q[2];
rz(-0.97270687) q[3];
sx q[3];
rz(-2.1199626) q[3];
sx q[3];
rz(-1.2620259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7369507) q[2];
sx q[2];
rz(-1.5664132) q[2];
sx q[2];
rz(2.289782) q[2];
rz(-2.8972054) q[3];
sx q[3];
rz(-2.4308725) q[3];
sx q[3];
rz(-1.3062668) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39563018) q[0];
sx q[0];
rz(-1.6599449) q[0];
sx q[0];
rz(-3.0356044) q[0];
rz(-1.7164187) q[1];
sx q[1];
rz(-1.9916078) q[1];
sx q[1];
rz(-1.0894159) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3857754) q[0];
sx q[0];
rz(-1.5013278) q[0];
sx q[0];
rz(-1.8088732) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0136497) q[2];
sx q[2];
rz(-2.2086655) q[2];
sx q[2];
rz(1.4952212) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8196667) q[1];
sx q[1];
rz(-1.2533256) q[1];
sx q[1];
rz(0.13369932) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.15481101) q[3];
sx q[3];
rz(-0.96199482) q[3];
sx q[3];
rz(1.7090635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.085422903) q[2];
sx q[2];
rz(-1.9499754) q[2];
sx q[2];
rz(2.5540111) q[2];
rz(-1.9116631) q[3];
sx q[3];
rz(-0.39837024) q[3];
sx q[3];
rz(2.5689094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9375482) q[0];
sx q[0];
rz(-1.9119104) q[0];
sx q[0];
rz(-2.6958534) q[0];
rz(2.2988689) q[1];
sx q[1];
rz(-0.80519599) q[1];
sx q[1];
rz(-2.3420948) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7460551) q[0];
sx q[0];
rz(-0.56395254) q[0];
sx q[0];
rz(2.3212295) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.92858814) q[2];
sx q[2];
rz(-2.1669496) q[2];
sx q[2];
rz(-0.10732574) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.362613) q[1];
sx q[1];
rz(-2.0019166) q[1];
sx q[1];
rz(-2.0126473) q[1];
x q[2];
rz(-1.4700674) q[3];
sx q[3];
rz(-0.58222773) q[3];
sx q[3];
rz(1.4034206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2004956) q[2];
sx q[2];
rz(-0.37797394) q[2];
sx q[2];
rz(-2.359158) q[2];
rz(-2.2109219) q[3];
sx q[3];
rz(-1.8698255) q[3];
sx q[3];
rz(-2.8455287) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3928669) q[0];
sx q[0];
rz(-0.94806945) q[0];
sx q[0];
rz(2.0178846) q[0];
rz(-2.3011082) q[1];
sx q[1];
rz(-1.1437462) q[1];
sx q[1];
rz(-1.3927654) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7531778) q[0];
sx q[0];
rz(-1.6196189) q[0];
sx q[0];
rz(-1.6724259) q[0];
rz(-pi) q[1];
x q[1];
rz(0.56864777) q[2];
sx q[2];
rz(-0.37867448) q[2];
sx q[2];
rz(-0.17998047) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9509134) q[1];
sx q[1];
rz(-1.2251405) q[1];
sx q[1];
rz(2.9596364) q[1];
x q[2];
rz(-0.61711611) q[3];
sx q[3];
rz(-1.6387306) q[3];
sx q[3];
rz(-1.9096494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.30354083) q[2];
sx q[2];
rz(-2.3736062) q[2];
sx q[2];
rz(1.6221907) q[2];
rz(2.0766025) q[3];
sx q[3];
rz(-1.0737373) q[3];
sx q[3];
rz(2.2860693) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72724718) q[0];
sx q[0];
rz(-2.1942744) q[0];
sx q[0];
rz(1.974768) q[0];
rz(0.30544063) q[1];
sx q[1];
rz(-2.0716397) q[1];
sx q[1];
rz(-2.6397612) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1239282) q[0];
sx q[0];
rz(-2.9249603) q[0];
sx q[0];
rz(2.6922035) q[0];
rz(-pi) q[1];
rz(2.3911693) q[2];
sx q[2];
rz(-2.1273251) q[2];
sx q[2];
rz(1.5897074) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8055685) q[1];
sx q[1];
rz(-2.7689287) q[1];
sx q[1];
rz(-1.3678846) q[1];
rz(1.0055528) q[3];
sx q[3];
rz(-2.0430419) q[3];
sx q[3];
rz(-0.22423553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.66494232) q[2];
sx q[2];
rz(-2.3822337) q[2];
sx q[2];
rz(-0.13258983) q[2];
rz(-2.8442123) q[3];
sx q[3];
rz(-1.5404276) q[3];
sx q[3];
rz(-0.69445777) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7731758) q[0];
sx q[0];
rz(-1.6650124) q[0];
sx q[0];
rz(1.52894) q[0];
rz(1.7932786) q[1];
sx q[1];
rz(-1.3765843) q[1];
sx q[1];
rz(-0.37337676) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0578846) q[0];
sx q[0];
rz(-2.7485131) q[0];
sx q[0];
rz(-1.1532093) q[0];
x q[1];
rz(-0.48252784) q[2];
sx q[2];
rz(-1.8410701) q[2];
sx q[2];
rz(-2.7729098) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8638525) q[1];
sx q[1];
rz(-1.2009632) q[1];
sx q[1];
rz(0.86994008) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4677603) q[3];
sx q[3];
rz(-2.8008409) q[3];
sx q[3];
rz(-2.1290523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1030964) q[2];
sx q[2];
rz(-0.080048397) q[2];
sx q[2];
rz(-0.45269629) q[2];
rz(-2.2321841) q[3];
sx q[3];
rz(-2.1287287) q[3];
sx q[3];
rz(-2.4694841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(1.7047003) q[0];
sx q[0];
rz(-1.1672651) q[0];
sx q[0];
rz(1.1080909) q[0];
rz(0.7946026) q[1];
sx q[1];
rz(-1.4789076) q[1];
sx q[1];
rz(2.2750003) q[1];
rz(1.5431326) q[2];
sx q[2];
rz(-2.4504708) q[2];
sx q[2];
rz(-2.8219555) q[2];
rz(-1.1232212) q[3];
sx q[3];
rz(-2.2991857) q[3];
sx q[3];
rz(-3.1274336) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
