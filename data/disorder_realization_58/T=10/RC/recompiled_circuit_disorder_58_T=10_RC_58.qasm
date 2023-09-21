OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9392202) q[0];
sx q[0];
rz(-0.4063172) q[0];
sx q[0];
rz(-2.321474) q[0];
rz(-0.36110538) q[1];
sx q[1];
rz(-2.5087924) q[1];
sx q[1];
rz(0.83067218) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3842073) q[0];
sx q[0];
rz(-1.877458) q[0];
sx q[0];
rz(2.7469809) q[0];
rz(-pi) q[1];
rz(-2.6248706) q[2];
sx q[2];
rz(-0.80846723) q[2];
sx q[2];
rz(0.39743039) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.1931397) q[1];
sx q[1];
rz(-1.8857191) q[1];
sx q[1];
rz(-0.83955168) q[1];
x q[2];
rz(-2.6942263) q[3];
sx q[3];
rz(-1.9379741) q[3];
sx q[3];
rz(-2.6255053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43710199) q[2];
sx q[2];
rz(-0.4318684) q[2];
sx q[2];
rz(0.1201771) q[2];
rz(1.9834571) q[3];
sx q[3];
rz(-1.7430256) q[3];
sx q[3];
rz(-0.91896287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08081089) q[0];
sx q[0];
rz(-1.7827543) q[0];
sx q[0];
rz(0.91180116) q[0];
rz(-2.3520825) q[1];
sx q[1];
rz(-0.98840886) q[1];
sx q[1];
rz(2.8149014) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17900285) q[0];
sx q[0];
rz(-2.0741182) q[0];
sx q[0];
rz(2.3110564) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9825282) q[2];
sx q[2];
rz(-2.0098364) q[2];
sx q[2];
rz(-0.95869267) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8596526) q[1];
sx q[1];
rz(-1.167206) q[1];
sx q[1];
rz(-2.2380026) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88300206) q[3];
sx q[3];
rz(-1.5044754) q[3];
sx q[3];
rz(-1.4955213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.5102753) q[2];
sx q[2];
rz(-2.3687506) q[2];
sx q[2];
rz(-1.2871683) q[2];
rz(3.0316947) q[3];
sx q[3];
rz(-1.7307614) q[3];
sx q[3];
rz(1.3818285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.592955) q[0];
sx q[0];
rz(-0.73919636) q[0];
sx q[0];
rz(2.8116995) q[0];
rz(-2.864481) q[1];
sx q[1];
rz(-1.3169293) q[1];
sx q[1];
rz(-2.0842016) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7085416) q[0];
sx q[0];
rz(-1.2305224) q[0];
sx q[0];
rz(0.021854594) q[0];
rz(1.7705275) q[2];
sx q[2];
rz(-1.235851) q[2];
sx q[2];
rz(-2.0470326) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.96506572) q[1];
sx q[1];
rz(-2.2245193) q[1];
sx q[1];
rz(2.7214126) q[1];
rz(1.3789165) q[3];
sx q[3];
rz(-1.9865611) q[3];
sx q[3];
rz(2.2172745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3588336) q[2];
sx q[2];
rz(-1.1153355) q[2];
sx q[2];
rz(1.3624297) q[2];
rz(2.5168915) q[3];
sx q[3];
rz(-2.0420859) q[3];
sx q[3];
rz(0.81956285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3574922) q[0];
sx q[0];
rz(-0.52629137) q[0];
sx q[0];
rz(-2.6065361) q[0];
rz(-1.1401945) q[1];
sx q[1];
rz(-1.3277206) q[1];
sx q[1];
rz(-0.16539703) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1918614) q[0];
sx q[0];
rz(-0.22338578) q[0];
sx q[0];
rz(0.97747691) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7521162) q[2];
sx q[2];
rz(-1.4099979) q[2];
sx q[2];
rz(-2.5374075) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0619229) q[1];
sx q[1];
rz(-1.8975782) q[1];
sx q[1];
rz(-2.5667739) q[1];
rz(-pi) q[2];
rz(1.9966647) q[3];
sx q[3];
rz(-0.59755675) q[3];
sx q[3];
rz(-1.715341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7455204) q[2];
sx q[2];
rz(-1.8053651) q[2];
sx q[2];
rz(-0.20544927) q[2];
rz(-2.0139587) q[3];
sx q[3];
rz(-1.1536359) q[3];
sx q[3];
rz(-1.2566465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66185343) q[0];
sx q[0];
rz(-0.91402188) q[0];
sx q[0];
rz(2.7868295) q[0];
rz(1.154249) q[1];
sx q[1];
rz(-0.92461363) q[1];
sx q[1];
rz(0.23194557) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1922798) q[0];
sx q[0];
rz(-0.67040196) q[0];
sx q[0];
rz(2.288726) q[0];
rz(-pi) q[1];
rz(-2.9734128) q[2];
sx q[2];
rz(-0.7025223) q[2];
sx q[2];
rz(2.0481734) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9310589) q[1];
sx q[1];
rz(-1.7909044) q[1];
sx q[1];
rz(-1.348043) q[1];
rz(-pi) q[2];
x q[2];
rz(0.36851818) q[3];
sx q[3];
rz(-1.5889865) q[3];
sx q[3];
rz(-2.8858678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0014687) q[2];
sx q[2];
rz(-1.8073558) q[2];
sx q[2];
rz(1.9011964) q[2];
rz(-2.5455348) q[3];
sx q[3];
rz(-1.3052992) q[3];
sx q[3];
rz(-1.3034472) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0444788) q[0];
sx q[0];
rz(-0.070274027) q[0];
sx q[0];
rz(-0.20275673) q[0];
rz(-0.98908201) q[1];
sx q[1];
rz(-1.443807) q[1];
sx q[1];
rz(0.99745497) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31842445) q[0];
sx q[0];
rz(-2.423375) q[0];
sx q[0];
rz(1.7757925) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1987655) q[2];
sx q[2];
rz(-0.87429201) q[2];
sx q[2];
rz(-0.55559413) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5207386) q[1];
sx q[1];
rz(-1.3315017) q[1];
sx q[1];
rz(-1.1362856) q[1];
x q[2];
rz(-2.2951179) q[3];
sx q[3];
rz(-1.9123565) q[3];
sx q[3];
rz(-0.70111707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9138907) q[2];
sx q[2];
rz(-1.9753549) q[2];
sx q[2];
rz(2.690199) q[2];
rz(2.732892) q[3];
sx q[3];
rz(-1.5390076) q[3];
sx q[3];
rz(1.9394978) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8354427) q[0];
sx q[0];
rz(-2.7375484) q[0];
sx q[0];
rz(0.62414449) q[0];
rz(1.6250601) q[1];
sx q[1];
rz(-0.25696483) q[1];
sx q[1];
rz(0.61378941) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81543024) q[0];
sx q[0];
rz(-2.3809732) q[0];
sx q[0];
rz(1.2994231) q[0];
x q[1];
rz(0.28366144) q[2];
sx q[2];
rz(-1.5093056) q[2];
sx q[2];
rz(-2.9279857) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1544513) q[1];
sx q[1];
rz(-0.11949355) q[1];
sx q[1];
rz(2.6277072) q[1];
rz(-0.11717637) q[3];
sx q[3];
rz(-1.1083318) q[3];
sx q[3];
rz(-1.6346491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7081786) q[2];
sx q[2];
rz(-0.91579473) q[2];
sx q[2];
rz(-1.9667352) q[2];
rz(-0.60837778) q[3];
sx q[3];
rz(-1.7374246) q[3];
sx q[3];
rz(1.7181989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90010086) q[0];
sx q[0];
rz(-1.2986203) q[0];
sx q[0];
rz(-0.4883782) q[0];
rz(1.6237367) q[1];
sx q[1];
rz(-1.7428215) q[1];
sx q[1];
rz(0.98446313) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0122089) q[0];
sx q[0];
rz(-2.1166271) q[0];
sx q[0];
rz(1.5270385) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.369957) q[2];
sx q[2];
rz(-1.2727591) q[2];
sx q[2];
rz(0.76921295) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9148548) q[1];
sx q[1];
rz(-1.9361155) q[1];
sx q[1];
rz(1.564333) q[1];
x q[2];
rz(-2.5823309) q[3];
sx q[3];
rz(-2.5945633) q[3];
sx q[3];
rz(0.83838851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.70665923) q[2];
sx q[2];
rz(-1.6985396) q[2];
sx q[2];
rz(2.6064176) q[2];
rz(2.0914071) q[3];
sx q[3];
rz(-1.8015367) q[3];
sx q[3];
rz(2.1323269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64131367) q[0];
sx q[0];
rz(-2.2135493) q[0];
sx q[0];
rz(0.6859268) q[0];
rz(-2.7507239) q[1];
sx q[1];
rz(-2.1978244) q[1];
sx q[1];
rz(-0.92591441) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029322421) q[0];
sx q[0];
rz(-1.3591213) q[0];
sx q[0];
rz(0.14365833) q[0];
x q[1];
rz(0.8530059) q[2];
sx q[2];
rz(-0.91663137) q[2];
sx q[2];
rz(-0.6086364) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1123062) q[1];
sx q[1];
rz(-0.46488133) q[1];
sx q[1];
rz(-1.9025365) q[1];
rz(-pi) q[2];
rz(1.4521452) q[3];
sx q[3];
rz(-1.3680397) q[3];
sx q[3];
rz(-2.6076536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.19568504) q[2];
sx q[2];
rz(-2.6699799) q[2];
sx q[2];
rz(-2.6055028) q[2];
rz(0.4195956) q[3];
sx q[3];
rz(-2.5511238) q[3];
sx q[3];
rz(1.6528116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0855899) q[0];
sx q[0];
rz(-0.35878006) q[0];
sx q[0];
rz(-2.7465903) q[0];
rz(1.5123873) q[1];
sx q[1];
rz(-1.2724266) q[1];
sx q[1];
rz(-2.1283456) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7618338) q[0];
sx q[0];
rz(-1.8031617) q[0];
sx q[0];
rz(-2.9985715) q[0];
rz(-3.0229438) q[2];
sx q[2];
rz(-0.35602202) q[2];
sx q[2];
rz(-2.7729386) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0777178) q[1];
sx q[1];
rz(-0.43242726) q[1];
sx q[1];
rz(0.23250154) q[1];
rz(-0.69443955) q[3];
sx q[3];
rz(-0.89530066) q[3];
sx q[3];
rz(-1.465786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.44832486) q[2];
sx q[2];
rz(-1.6222745) q[2];
sx q[2];
rz(0.75941336) q[2];
rz(1.365472) q[3];
sx q[3];
rz(-1.1080192) q[3];
sx q[3];
rz(0.56994462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7286745) q[0];
sx q[0];
rz(-1.2644132) q[0];
sx q[0];
rz(1.2046474) q[0];
rz(0.57327523) q[1];
sx q[1];
rz(-1.9075867) q[1];
sx q[1];
rz(1.8214068) q[1];
rz(-1.7309932) q[2];
sx q[2];
rz(-2.6916531) q[2];
sx q[2];
rz(2.0315363) q[2];
rz(-1.4952954) q[3];
sx q[3];
rz(-2.4423238) q[3];
sx q[3];
rz(-2.8764976) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];