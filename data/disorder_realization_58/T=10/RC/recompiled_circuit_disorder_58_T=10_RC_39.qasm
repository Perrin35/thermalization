OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.20237246) q[0];
sx q[0];
rz(-2.7352754) q[0];
sx q[0];
rz(2.321474) q[0];
rz(-0.36110538) q[1];
sx q[1];
rz(0.63280025) q[1];
sx q[1];
rz(11.735698) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3842073) q[0];
sx q[0];
rz(-1.877458) q[0];
sx q[0];
rz(-0.39461179) q[0];
rz(-pi) q[1];
rz(2.0482424) q[2];
sx q[2];
rz(-0.89077836) q[2];
sx q[2];
rz(-2.8505461) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1066061) q[1];
sx q[1];
rz(-2.258746) q[1];
sx q[1];
rz(-0.41253849) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6942263) q[3];
sx q[3];
rz(-1.2036185) q[3];
sx q[3];
rz(-2.6255053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.43710199) q[2];
sx q[2];
rz(-0.4318684) q[2];
sx q[2];
rz(-0.1201771) q[2];
rz(-1.1581356) q[3];
sx q[3];
rz(-1.3985671) q[3];
sx q[3];
rz(0.91896287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0607818) q[0];
sx q[0];
rz(-1.7827543) q[0];
sx q[0];
rz(-0.91180116) q[0];
rz(0.78951019) q[1];
sx q[1];
rz(-2.1531838) q[1];
sx q[1];
rz(-2.8149014) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3347496) q[0];
sx q[0];
rz(-2.202889) q[0];
sx q[0];
rz(2.500781) q[0];
rz(-pi) q[1];
rz(-2.9825282) q[2];
sx q[2];
rz(-2.0098364) q[2];
sx q[2];
rz(-0.95869267) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3152299) q[1];
sx q[1];
rz(-0.76347199) q[1];
sx q[1];
rz(0.96674322) q[1];
x q[2];
rz(1.4665524) q[3];
sx q[3];
rz(-0.6904656) q[3];
sx q[3];
rz(-0.0052099293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6313173) q[2];
sx q[2];
rz(-0.77284208) q[2];
sx q[2];
rz(-1.8544244) q[2];
rz(-0.10989799) q[3];
sx q[3];
rz(-1.7307614) q[3];
sx q[3];
rz(-1.7597642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54863769) q[0];
sx q[0];
rz(-2.4023963) q[0];
sx q[0];
rz(-2.8116995) q[0];
rz(2.864481) q[1];
sx q[1];
rz(-1.3169293) q[1];
sx q[1];
rz(2.0842016) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7085416) q[0];
sx q[0];
rz(-1.2305224) q[0];
sx q[0];
rz(-3.1197381) q[0];
rz(-0.51809394) q[2];
sx q[2];
rz(-2.7535541) q[2];
sx q[2];
rz(-0.54259091) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.33144618) q[1];
sx q[1];
rz(-2.3815037) q[1];
sx q[1];
rz(2.0600832) q[1];
rz(-pi) q[2];
rz(-2.7338545) q[3];
sx q[3];
rz(-2.6860415) q[3];
sx q[3];
rz(1.7689442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.7827591) q[2];
sx q[2];
rz(-1.1153355) q[2];
sx q[2];
rz(1.7791629) q[2];
rz(2.5168915) q[3];
sx q[3];
rz(-1.0995068) q[3];
sx q[3];
rz(-0.81956285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3574922) q[0];
sx q[0];
rz(-0.52629137) q[0];
sx q[0];
rz(0.5350565) q[0];
rz(-1.1401945) q[1];
sx q[1];
rz(-1.813872) q[1];
sx q[1];
rz(0.16539703) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3446942) q[0];
sx q[0];
rz(-1.3860774) q[0];
sx q[0];
rz(3.0152507) q[0];
x q[1];
rz(-2.9781614) q[2];
sx q[2];
rz(-1.7497517) q[2];
sx q[2];
rz(-2.2043259) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.71388984) q[1];
sx q[1];
rz(-2.1117003) q[1];
sx q[1];
rz(-1.9546024) q[1];
rz(1.0159675) q[3];
sx q[3];
rz(-1.3362243) q[3];
sx q[3];
rz(0.50333422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.39607221) q[2];
sx q[2];
rz(-1.3362276) q[2];
sx q[2];
rz(0.20544927) q[2];
rz(2.0139587) q[3];
sx q[3];
rz(-1.9879568) q[3];
sx q[3];
rz(1.8849461) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66185343) q[0];
sx q[0];
rz(-2.2275708) q[0];
sx q[0];
rz(2.7868295) q[0];
rz(-1.9873437) q[1];
sx q[1];
rz(-2.216979) q[1];
sx q[1];
rz(2.9096471) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97868279) q[0];
sx q[0];
rz(-1.1497578) q[0];
sx q[0];
rz(1.0324423) q[0];
rz(-pi) q[1];
rz(-1.7115713) q[2];
sx q[2];
rz(-0.88015926) q[2];
sx q[2];
rz(1.8292793) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2105337) q[1];
sx q[1];
rz(-1.7909044) q[1];
sx q[1];
rz(1.348043) q[1];
rz(-0.050458126) q[3];
sx q[3];
rz(-0.36894635) q[3];
sx q[3];
rz(1.7794533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.140124) q[2];
sx q[2];
rz(-1.8073558) q[2];
sx q[2];
rz(1.2403963) q[2];
rz(-2.5455348) q[3];
sx q[3];
rz(-1.8362935) q[3];
sx q[3];
rz(1.3034472) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0444788) q[0];
sx q[0];
rz(-3.0713186) q[0];
sx q[0];
rz(0.20275673) q[0];
rz(-2.1525106) q[1];
sx q[1];
rz(-1.6977856) q[1];
sx q[1];
rz(0.99745497) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8231682) q[0];
sx q[0];
rz(-2.423375) q[0];
sx q[0];
rz(1.3658001) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41001292) q[2];
sx q[2];
rz(-2.3668681) q[2];
sx q[2];
rz(2.0395525) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5207386) q[1];
sx q[1];
rz(-1.8100909) q[1];
sx q[1];
rz(1.1362856) q[1];
rz(2.698425) q[3];
sx q[3];
rz(-0.89649761) q[3];
sx q[3];
rz(0.5815732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.22770195) q[2];
sx q[2];
rz(-1.9753549) q[2];
sx q[2];
rz(-2.690199) q[2];
rz(-0.40870062) q[3];
sx q[3];
rz(-1.6025851) q[3];
sx q[3];
rz(1.2020948) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8354427) q[0];
sx q[0];
rz(-0.4040443) q[0];
sx q[0];
rz(0.62414449) q[0];
rz(-1.5165326) q[1];
sx q[1];
rz(-0.25696483) q[1];
sx q[1];
rz(-2.5278032) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44871556) q[0];
sx q[0];
rz(-2.2971417) q[0];
sx q[0];
rz(-2.8918299) q[0];
x q[1];
rz(-1.6348398) q[2];
sx q[2];
rz(-1.8539068) q[2];
sx q[2];
rz(-1.7664906) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.63749667) q[1];
sx q[1];
rz(-1.4667965) q[1];
sx q[1];
rz(-1.6297479) q[1];
rz(-pi) q[2];
rz(1.1055787) q[3];
sx q[3];
rz(-1.675616) q[3];
sx q[3];
rz(-0.011381586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7081786) q[2];
sx q[2];
rz(-0.91579473) q[2];
sx q[2];
rz(1.9667352) q[2];
rz(2.5332149) q[3];
sx q[3];
rz(-1.404168) q[3];
sx q[3];
rz(-1.7181989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2414918) q[0];
sx q[0];
rz(-1.2986203) q[0];
sx q[0];
rz(-0.4883782) q[0];
rz(-1.6237367) q[1];
sx q[1];
rz(-1.7428215) q[1];
sx q[1];
rz(2.1571295) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1293837) q[0];
sx q[0];
rz(-2.1166271) q[0];
sx q[0];
rz(1.6145541) q[0];
rz(1.1659053) q[2];
sx q[2];
rz(-2.3003909) q[2];
sx q[2];
rz(1.0798432) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8967646) q[1];
sx q[1];
rz(-0.36537376) q[1];
sx q[1];
rz(0.016896292) q[1];
rz(0.55926178) q[3];
sx q[3];
rz(-0.54702938) q[3];
sx q[3];
rz(2.3032041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4349334) q[2];
sx q[2];
rz(-1.443053) q[2];
sx q[2];
rz(2.6064176) q[2];
rz(-2.0914071) q[3];
sx q[3];
rz(-1.340056) q[3];
sx q[3];
rz(-1.0092658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.500279) q[0];
sx q[0];
rz(-0.92804337) q[0];
sx q[0];
rz(0.6859268) q[0];
rz(0.39086875) q[1];
sx q[1];
rz(-2.1978244) q[1];
sx q[1];
rz(-0.92591441) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5092963) q[0];
sx q[0];
rz(-2.8863781) q[0];
sx q[0];
rz(2.1584828) q[0];
rz(2.2885867) q[2];
sx q[2];
rz(-0.91663137) q[2];
sx q[2];
rz(0.6086364) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8402108) q[1];
sx q[1];
rz(-1.4242607) q[1];
sx q[1];
rz(-2.0135897) q[1];
rz(1.4521452) q[3];
sx q[3];
rz(-1.7735529) q[3];
sx q[3];
rz(-0.53393902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9459076) q[2];
sx q[2];
rz(-2.6699799) q[2];
sx q[2];
rz(-2.6055028) q[2];
rz(-0.4195956) q[3];
sx q[3];
rz(-0.59046888) q[3];
sx q[3];
rz(1.6528116) q[3];
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
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0560028) q[0];
sx q[0];
rz(-2.7828126) q[0];
sx q[0];
rz(-2.7465903) q[0];
rz(1.5123873) q[1];
sx q[1];
rz(-1.2724266) q[1];
sx q[1];
rz(1.013247) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2241867) q[0];
sx q[0];
rz(-1.7099483) q[0];
sx q[0];
rz(-1.805472) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5268065) q[2];
sx q[2];
rz(-1.9242052) q[2];
sx q[2];
rz(-2.8994438) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0638748) q[1];
sx q[1];
rz(-2.7091654) q[1];
sx q[1];
rz(-2.9090911) q[1];
rz(-2.3771044) q[3];
sx q[3];
rz(-1.0478684) q[3];
sx q[3];
rz(-0.58512277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.44832486) q[2];
sx q[2];
rz(-1.6222745) q[2];
sx q[2];
rz(2.3821793) q[2];
rz(-1.365472) q[3];
sx q[3];
rz(-2.0335734) q[3];
sx q[3];
rz(0.56994462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7286745) q[0];
sx q[0];
rz(-1.8771794) q[0];
sx q[0];
rz(-1.9369453) q[0];
rz(-2.5683174) q[1];
sx q[1];
rz(-1.9075867) q[1];
sx q[1];
rz(1.8214068) q[1];
rz(1.7309932) q[2];
sx q[2];
rz(-0.44993958) q[2];
sx q[2];
rz(-1.1100563) q[2];
rz(3.0782386) q[3];
sx q[3];
rz(-2.2676716) q[3];
sx q[3];
rz(0.36361658) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];