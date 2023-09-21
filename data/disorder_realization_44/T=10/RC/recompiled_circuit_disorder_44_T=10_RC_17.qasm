OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2919579) q[0];
sx q[0];
rz(-2.7014974) q[0];
sx q[0];
rz(-0.13719288) q[0];
rz(1.4057012) q[1];
sx q[1];
rz(4.5448137) q[1];
sx q[1];
rz(9.9546976) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0041381) q[0];
sx q[0];
rz(-1.9507427) q[0];
sx q[0];
rz(0.11560346) q[0];
x q[1];
rz(-2.797384) q[2];
sx q[2];
rz(-1.95103) q[2];
sx q[2];
rz(0.76438475) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4094028) q[1];
sx q[1];
rz(-2.4362872) q[1];
sx q[1];
rz(0.7440872) q[1];
x q[2];
rz(2.3832541) q[3];
sx q[3];
rz(-1.3876649) q[3];
sx q[3];
rz(1.5116215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.68937504) q[2];
sx q[2];
rz(-1.3000501) q[2];
sx q[2];
rz(-0.33660647) q[2];
rz(-1.6254788) q[3];
sx q[3];
rz(-0.55364004) q[3];
sx q[3];
rz(1.6158993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34823725) q[0];
sx q[0];
rz(-2.0331148) q[0];
sx q[0];
rz(3.120378) q[0];
rz(-1.9477828) q[1];
sx q[1];
rz(-2.1021011) q[1];
sx q[1];
rz(-0.83591998) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6244038) q[0];
sx q[0];
rz(-1.7096585) q[0];
sx q[0];
rz(-0.018272321) q[0];
rz(-pi) q[1];
rz(0.40936562) q[2];
sx q[2];
rz(-1.0592807) q[2];
sx q[2];
rz(-0.53586938) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1509717) q[1];
sx q[1];
rz(-2.4358106) q[1];
sx q[1];
rz(-1.1433931) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1777982) q[3];
sx q[3];
rz(-0.8702232) q[3];
sx q[3];
rz(-2.8515479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8643643) q[2];
sx q[2];
rz(-1.1479062) q[2];
sx q[2];
rz(1.7956087) q[2];
rz(2.7820382) q[3];
sx q[3];
rz(-0.94272009) q[3];
sx q[3];
rz(2.6446222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42831746) q[0];
sx q[0];
rz(-2.064216) q[0];
sx q[0];
rz(-1.0536449) q[0];
rz(1.9127649) q[1];
sx q[1];
rz(-1.6002974) q[1];
sx q[1];
rz(2.704481) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8374098) q[0];
sx q[0];
rz(-1.8619616) q[0];
sx q[0];
rz(3.0512179) q[0];
x q[1];
rz(1.978546) q[2];
sx q[2];
rz(-1.7951269) q[2];
sx q[2];
rz(-0.11220223) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5324085) q[1];
sx q[1];
rz(-0.44645616) q[1];
sx q[1];
rz(2.8162454) q[1];
rz(1.6821074) q[3];
sx q[3];
rz(-2.6442332) q[3];
sx q[3];
rz(1.5670083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.019471021) q[2];
sx q[2];
rz(-0.78142587) q[2];
sx q[2];
rz(-1.0220698) q[2];
rz(-1.2381037) q[3];
sx q[3];
rz(-2.759203) q[3];
sx q[3];
rz(2.7220272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6195174) q[0];
sx q[0];
rz(-1.8958805) q[0];
sx q[0];
rz(-2.1602901) q[0];
rz(-0.13521067) q[1];
sx q[1];
rz(-2.0573261) q[1];
sx q[1];
rz(-0.19128004) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0606196) q[0];
sx q[0];
rz(-1.7195716) q[0];
sx q[0];
rz(3.0540375) q[0];
rz(1.26881) q[2];
sx q[2];
rz(-0.75690818) q[2];
sx q[2];
rz(0.41727558) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8297255) q[1];
sx q[1];
rz(-0.79677454) q[1];
sx q[1];
rz(1.4273248) q[1];
rz(2.6077765) q[3];
sx q[3];
rz(-1.2663519) q[3];
sx q[3];
rz(-0.75418762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.68025756) q[2];
sx q[2];
rz(-2.1562083) q[2];
sx q[2];
rz(-2.130924) q[2];
rz(0.7615532) q[3];
sx q[3];
rz(-1.9617617) q[3];
sx q[3];
rz(-0.23553577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.8028832) q[0];
sx q[0];
rz(-0.25512472) q[0];
sx q[0];
rz(-0.55661911) q[0];
rz(-3.026475) q[1];
sx q[1];
rz(-1.8042253) q[1];
sx q[1];
rz(-0.97250485) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8079677) q[0];
sx q[0];
rz(-1.4687612) q[0];
sx q[0];
rz(-1.6385965) q[0];
rz(0.33072492) q[2];
sx q[2];
rz(-0.84825584) q[2];
sx q[2];
rz(-1.6387788) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3481808) q[1];
sx q[1];
rz(-0.23985292) q[1];
sx q[1];
rz(-0.94437771) q[1];
x q[2];
rz(0.036690849) q[3];
sx q[3];
rz(-0.9551691) q[3];
sx q[3];
rz(-0.43945593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3187023) q[2];
sx q[2];
rz(-1.0914785) q[2];
sx q[2];
rz(1.548432) q[2];
rz(-1.3657773) q[3];
sx q[3];
rz(-2.8184991) q[3];
sx q[3];
rz(0.95388609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71972972) q[0];
sx q[0];
rz(-1.2635764) q[0];
sx q[0];
rz(1.7156037) q[0];
rz(-1.0643719) q[1];
sx q[1];
rz(-2.1247037) q[1];
sx q[1];
rz(-0.37429601) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4842589) q[0];
sx q[0];
rz(-1.1614292) q[0];
sx q[0];
rz(-1.2249468) q[0];
rz(-1.3099953) q[2];
sx q[2];
rz(-1.0909832) q[2];
sx q[2];
rz(0.11670437) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1230871) q[1];
sx q[1];
rz(-1.0235041) q[1];
sx q[1];
rz(1.2612543) q[1];
x q[2];
rz(-2.4738594) q[3];
sx q[3];
rz(-0.99620921) q[3];
sx q[3];
rz(-2.7001065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8950243) q[2];
sx q[2];
rz(-0.66528577) q[2];
sx q[2];
rz(-2.1833615) q[2];
rz(0.22917497) q[3];
sx q[3];
rz(-1.6848247) q[3];
sx q[3];
rz(2.5206101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7763057) q[0];
sx q[0];
rz(-1.1927274) q[0];
sx q[0];
rz(-0.90674415) q[0];
rz(-2.0523741) q[1];
sx q[1];
rz(-1.6420495) q[1];
sx q[1];
rz(1.3100756) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7999254) q[0];
sx q[0];
rz(-2.7503715) q[0];
sx q[0];
rz(-0.13473405) q[0];
rz(-pi) q[1];
rz(-1.9016978) q[2];
sx q[2];
rz(-0.885193) q[2];
sx q[2];
rz(2.9943525) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.20565614) q[1];
sx q[1];
rz(-1.2381136) q[1];
sx q[1];
rz(1.3693621) q[1];
rz(-pi) q[2];
rz(1.5433611) q[3];
sx q[3];
rz(-2.1260288) q[3];
sx q[3];
rz(2.1813986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.92581302) q[2];
sx q[2];
rz(-0.44162193) q[2];
sx q[2];
rz(-1.4833935) q[2];
rz(0.27967927) q[3];
sx q[3];
rz(-2.1488583) q[3];
sx q[3];
rz(3.083995) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9782372) q[0];
sx q[0];
rz(-1.6219448) q[0];
sx q[0];
rz(-0.21959198) q[0];
rz(0.50312463) q[1];
sx q[1];
rz(-2.2527835) q[1];
sx q[1];
rz(0.84987744) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4445164) q[0];
sx q[0];
rz(-2.8651617) q[0];
sx q[0];
rz(-0.96093775) q[0];
rz(-pi) q[1];
rz(-1.4555898) q[2];
sx q[2];
rz(-2.4552279) q[2];
sx q[2];
rz(2.3442868) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1666607) q[1];
sx q[1];
rz(-0.76172511) q[1];
sx q[1];
rz(1.6208266) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0073118) q[3];
sx q[3];
rz(-2.8562162) q[3];
sx q[3];
rz(0.34261045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.59051096) q[2];
sx q[2];
rz(-1.8293646) q[2];
sx q[2];
rz(-1.760651) q[2];
rz(-0.75602174) q[3];
sx q[3];
rz(-2.9383926) q[3];
sx q[3];
rz(0.35593629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6324156) q[0];
sx q[0];
rz(-2.248705) q[0];
sx q[0];
rz(2.7365622) q[0];
rz(-0.45267725) q[1];
sx q[1];
rz(-2.15937) q[1];
sx q[1];
rz(-1.2776432) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15690878) q[0];
sx q[0];
rz(-0.55159969) q[0];
sx q[0];
rz(-2.6916235) q[0];
rz(-pi) q[1];
rz(0.69127609) q[2];
sx q[2];
rz(-1.9553767) q[2];
sx q[2];
rz(-0.28011766) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9459878) q[1];
sx q[1];
rz(-1.1254278) q[1];
sx q[1];
rz(-3.1006378) q[1];
rz(-pi) q[2];
rz(1.9926662) q[3];
sx q[3];
rz(-2.359169) q[3];
sx q[3];
rz(-2.9187834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3383125) q[2];
sx q[2];
rz(-2.7351604) q[2];
sx q[2];
rz(-2.7837616) q[2];
rz(1.7221649) q[3];
sx q[3];
rz(-1.2737041) q[3];
sx q[3];
rz(-1.0740124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91838592) q[0];
sx q[0];
rz(-0.077843852) q[0];
sx q[0];
rz(-3.0293368) q[0];
rz(-0.90011251) q[1];
sx q[1];
rz(-2.0745514) q[1];
sx q[1];
rz(2.9311438) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68073273) q[0];
sx q[0];
rz(-1.0675149) q[0];
sx q[0];
rz(1.7778648) q[0];
rz(-2.9774882) q[2];
sx q[2];
rz(-1.9774984) q[2];
sx q[2];
rz(-1.1526398) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3810972) q[1];
sx q[1];
rz(-1.3638745) q[1];
sx q[1];
rz(-0.31462545) q[1];
x q[2];
rz(-2.3841303) q[3];
sx q[3];
rz(-2.8031581) q[3];
sx q[3];
rz(-0.38871845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.18008733) q[2];
sx q[2];
rz(-0.6663565) q[2];
sx q[2];
rz(1.5562742) q[2];
rz(-1.8680343) q[3];
sx q[3];
rz(-0.62265101) q[3];
sx q[3];
rz(0.20726985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56959854) q[0];
sx q[0];
rz(-0.8710237) q[0];
sx q[0];
rz(-1.3652753) q[0];
rz(0.81644425) q[1];
sx q[1];
rz(-1.888231) q[1];
sx q[1];
rz(2.9838557) q[1];
rz(2.9755637) q[2];
sx q[2];
rz(-1.1062853) q[2];
sx q[2];
rz(1.7044978) q[2];
rz(2.7995085) q[3];
sx q[3];
rz(-2.2707006) q[3];
sx q[3];
rz(-2.0778098) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];