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
rz(-2.8357808) q[0];
sx q[0];
rz(-0.6337136) q[0];
sx q[0];
rz(-0.60854882) q[0];
rz(2.4701056) q[1];
sx q[1];
rz(-0.66317135) q[1];
sx q[1];
rz(-1.4594249) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08069399) q[0];
sx q[0];
rz(-2.0804724) q[0];
sx q[0];
rz(-2.8721548) q[0];
rz(2.0107277) q[2];
sx q[2];
rz(-2.6509883) q[2];
sx q[2];
rz(-3.0705796) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1009095) q[1];
sx q[1];
rz(-0.95129943) q[1];
sx q[1];
rz(-0.49609523) q[1];
rz(-pi) q[2];
rz(0.75992775) q[3];
sx q[3];
rz(-1.4455405) q[3];
sx q[3];
rz(-1.6881936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3964316) q[2];
sx q[2];
rz(-0.67239422) q[2];
sx q[2];
rz(2.3154837) q[2];
rz(-0.36469001) q[3];
sx q[3];
rz(-1.5267905) q[3];
sx q[3];
rz(2.8472692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.9593338) q[0];
sx q[0];
rz(-1.8987645) q[0];
sx q[0];
rz(2.6958595) q[0];
rz(0.590473) q[1];
sx q[1];
rz(-1.0346552) q[1];
sx q[1];
rz(0.42764923) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2954461) q[0];
sx q[0];
rz(-1.5238191) q[0];
sx q[0];
rz(2.0278005) q[0];
rz(-2.325133) q[2];
sx q[2];
rz(-1.1414981) q[2];
sx q[2];
rz(0.72390899) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.44438293) q[1];
sx q[1];
rz(-2.1661609) q[1];
sx q[1];
rz(-0.31636379) q[1];
rz(-pi) q[2];
rz(-0.0049875445) q[3];
sx q[3];
rz(-0.79646275) q[3];
sx q[3];
rz(2.4465268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6815971) q[2];
sx q[2];
rz(-2.5619016) q[2];
sx q[2];
rz(1.0350636) q[2];
rz(1.5929818) q[3];
sx q[3];
rz(-0.17954738) q[3];
sx q[3];
rz(-0.88919324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11397938) q[0];
sx q[0];
rz(-3.1359105) q[0];
sx q[0];
rz(-0.11153829) q[0];
rz(-1.5533718) q[1];
sx q[1];
rz(-0.40451834) q[1];
sx q[1];
rz(-2.9389971) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83765471) q[0];
sx q[0];
rz(-1.1753774) q[0];
sx q[0];
rz(2.9962792) q[0];
x q[1];
rz(1.2738704) q[2];
sx q[2];
rz(-1.8477173) q[2];
sx q[2];
rz(0.94147791) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2592273) q[1];
sx q[1];
rz(-2.276032) q[1];
sx q[1];
rz(2.3036876) q[1];
rz(-0.58835619) q[3];
sx q[3];
rz(-0.6783456) q[3];
sx q[3];
rz(2.7916186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.69278875) q[2];
sx q[2];
rz(-2.0710129) q[2];
sx q[2];
rz(-2.1602574) q[2];
rz(-0.28686178) q[3];
sx q[3];
rz(-2.562264) q[3];
sx q[3];
rz(-2.9441492) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9051179) q[0];
sx q[0];
rz(-0.050967) q[0];
sx q[0];
rz(0.68848759) q[0];
rz(-0.14432898) q[1];
sx q[1];
rz(-1.9537484) q[1];
sx q[1];
rz(-2.3932638) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8546974) q[0];
sx q[0];
rz(-2.4208491) q[0];
sx q[0];
rz(0.067528226) q[0];
rz(-pi) q[1];
rz(0.26864895) q[2];
sx q[2];
rz(-1.2792818) q[2];
sx q[2];
rz(-0.47713258) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6627745) q[1];
sx q[1];
rz(-1.2726901) q[1];
sx q[1];
rz(-1.3762022) q[1];
rz(1.9311936) q[3];
sx q[3];
rz(-2.4080417) q[3];
sx q[3];
rz(-1.4195132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0904842) q[2];
sx q[2];
rz(-1.8277617) q[2];
sx q[2];
rz(-1.8947821) q[2];
rz(3.0221853) q[3];
sx q[3];
rz(-1.2179008) q[3];
sx q[3];
rz(-2.8363805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3731821) q[0];
sx q[0];
rz(-2.7550582) q[0];
sx q[0];
rz(-0.61491948) q[0];
rz(2.3070295) q[1];
sx q[1];
rz(-1.395023) q[1];
sx q[1];
rz(0.83647299) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1008719) q[0];
sx q[0];
rz(-3.1060069) q[0];
sx q[0];
rz(-1.3077223) q[0];
rz(-pi) q[1];
rz(2.1600989) q[2];
sx q[2];
rz(-0.98926614) q[2];
sx q[2];
rz(-2.6133693) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.78042049) q[1];
sx q[1];
rz(-0.57406146) q[1];
sx q[1];
rz(-1.7758382) q[1];
x q[2];
rz(-2.3616124) q[3];
sx q[3];
rz(-2.2298919) q[3];
sx q[3];
rz(2.7798228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5323083) q[2];
sx q[2];
rz(-0.60695761) q[2];
sx q[2];
rz(2.2212846) q[2];
rz(2.3406384) q[3];
sx q[3];
rz(-2.6194173) q[3];
sx q[3];
rz(-2.8029158) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1717218) q[0];
sx q[0];
rz(-2.0750177) q[0];
sx q[0];
rz(-0.26611662) q[0];
rz(-1.7099821) q[1];
sx q[1];
rz(-2.0224729) q[1];
sx q[1];
rz(-2.5439579) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1640329) q[0];
sx q[0];
rz(-2.4983239) q[0];
sx q[0];
rz(-2.8215762) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.581988) q[2];
sx q[2];
rz(-1.3897685) q[2];
sx q[2];
rz(-2.343296) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8460076) q[1];
sx q[1];
rz(-0.69900988) q[1];
sx q[1];
rz(-3/(11*pi)) q[1];
rz(-1.131874) q[3];
sx q[3];
rz(-2.7898129) q[3];
sx q[3];
rz(-1.9466512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0133682) q[2];
sx q[2];
rz(-2.8743663) q[2];
sx q[2];
rz(2.5426148) q[2];
rz(-2.525575) q[3];
sx q[3];
rz(-0.44121656) q[3];
sx q[3];
rz(1.0429355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1004341) q[0];
sx q[0];
rz(-2.2293595) q[0];
sx q[0];
rz(2.8225733) q[0];
rz(3.0306385) q[1];
sx q[1];
rz(-1.3919421) q[1];
sx q[1];
rz(-0.14161938) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9931147) q[0];
sx q[0];
rz(-1.6671868) q[0];
sx q[0];
rz(-0.68401159) q[0];
x q[1];
rz(-1.0194112) q[2];
sx q[2];
rz(-1.5670318) q[2];
sx q[2];
rz(0.55094811) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7748877) q[1];
sx q[1];
rz(-2.1122518) q[1];
sx q[1];
rz(-1.790994) q[1];
rz(-pi) q[2];
rz(1.3124489) q[3];
sx q[3];
rz(-0.53505361) q[3];
sx q[3];
rz(2.759397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5854599) q[2];
sx q[2];
rz(-1.0254878) q[2];
sx q[2];
rz(0.11036135) q[2];
rz(-0.50802556) q[3];
sx q[3];
rz(-0.042162687) q[3];
sx q[3];
rz(1.149811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13462774) q[0];
sx q[0];
rz(-2.4899794) q[0];
sx q[0];
rz(1.2192669) q[0];
rz(-0.02267516) q[1];
sx q[1];
rz(-2.2494648) q[1];
sx q[1];
rz(2.5882744) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5641877) q[0];
sx q[0];
rz(-2.4779442) q[0];
sx q[0];
rz(1.7148561) q[0];
x q[1];
rz(-0.99626096) q[2];
sx q[2];
rz(-1.8991514) q[2];
sx q[2];
rz(2.7789214) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.07456116) q[1];
sx q[1];
rz(-2.0843049) q[1];
sx q[1];
rz(1.1719955) q[1];
x q[2];
rz(1.2723921) q[3];
sx q[3];
rz(-1.8074028) q[3];
sx q[3];
rz(-2.2158156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.40552178) q[2];
sx q[2];
rz(-2.5260479) q[2];
sx q[2];
rz(1.3496529) q[2];
rz(0.060529709) q[3];
sx q[3];
rz(-0.90977257) q[3];
sx q[3];
rz(-0.6282261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
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
rz(-2.8796006) q[0];
sx q[0];
rz(-2.1475726) q[0];
sx q[0];
rz(-3.1349728) q[0];
rz(0.63893843) q[1];
sx q[1];
rz(-2.9220351) q[1];
sx q[1];
rz(-2.7117597) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2868067) q[0];
sx q[0];
rz(-1.2116951) q[0];
sx q[0];
rz(0.39639985) q[0];
x q[1];
rz(2.8877385) q[2];
sx q[2];
rz(-1.1510669) q[2];
sx q[2];
rz(-1.8373289) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.92569578) q[1];
sx q[1];
rz(-1.5491423) q[1];
sx q[1];
rz(-2.3587061) q[1];
rz(-pi) q[2];
rz(2.5561629) q[3];
sx q[3];
rz(-1.5133891) q[3];
sx q[3];
rz(1.9229744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2868353) q[2];
sx q[2];
rz(-1.9660051) q[2];
sx q[2];
rz(0.01469928) q[2];
rz(-2.0059351) q[3];
sx q[3];
rz(-1.1123927) q[3];
sx q[3];
rz(1.3057825) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.056674615) q[0];
sx q[0];
rz(-0.25923964) q[0];
sx q[0];
rz(1.8490476) q[0];
rz(-2.9855285) q[1];
sx q[1];
rz(-0.81764692) q[1];
sx q[1];
rz(2.3110716) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9915402) q[0];
sx q[0];
rz(-2.4084598) q[0];
sx q[0];
rz(-1.7197648) q[0];
rz(1.8761329) q[2];
sx q[2];
rz(-2.3407901) q[2];
sx q[2];
rz(2.3126471) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3519684) q[1];
sx q[1];
rz(-1.3543173) q[1];
sx q[1];
rz(-0.08190544) q[1];
rz(-pi) q[2];
x q[2];
rz(0.36673185) q[3];
sx q[3];
rz(-1.368513) q[3];
sx q[3];
rz(-2.3679581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0326651) q[2];
sx q[2];
rz(-0.69643164) q[2];
sx q[2];
rz(2.6984974) q[2];
rz(-2.9923934) q[3];
sx q[3];
rz(-2.7982893) q[3];
sx q[3];
rz(-2.8086015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14015848) q[0];
sx q[0];
rz(-0.59068155) q[0];
sx q[0];
rz(-0.97472192) q[0];
rz(-2.1433266) q[1];
sx q[1];
rz(-1.9010192) q[1];
sx q[1];
rz(-0.98767282) q[1];
rz(-1.3692296) q[2];
sx q[2];
rz(-1.907363) q[2];
sx q[2];
rz(3.0628352) q[2];
rz(-0.90726388) q[3];
sx q[3];
rz(-1.2259007) q[3];
sx q[3];
rz(-2.905734) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
