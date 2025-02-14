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
rz(2.7414275) q[0];
sx q[0];
rz(-2.6097809) q[0];
sx q[0];
rz(-2.3398633) q[0];
rz(0.62893686) q[1];
sx q[1];
rz(-1.8328272) q[1];
sx q[1];
rz(-2.7523249) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9297732) q[0];
sx q[0];
rz(-2.4113829) q[0];
sx q[0];
rz(0.41359652) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9603329) q[2];
sx q[2];
rz(-1.6158293) q[2];
sx q[2];
rz(-0.35895106) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9323665) q[1];
sx q[1];
rz(-2.3079909) q[1];
sx q[1];
rz(-0.7070138) q[1];
x q[2];
rz(0.84343976) q[3];
sx q[3];
rz(-1.5411845) q[3];
sx q[3];
rz(0.51472291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6486711) q[2];
sx q[2];
rz(-0.45009437) q[2];
sx q[2];
rz(-2.8802803) q[2];
rz(0.35604769) q[3];
sx q[3];
rz(-1.1084403) q[3];
sx q[3];
rz(1.0935008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4556731) q[0];
sx q[0];
rz(-2.5643667) q[0];
sx q[0];
rz(2.8976231) q[0];
rz(0.4295373) q[1];
sx q[1];
rz(-1.669084) q[1];
sx q[1];
rz(-0.74091774) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8888433) q[0];
sx q[0];
rz(-2.6359632) q[0];
sx q[0];
rz(0.83871418) q[0];
rz(1.6119611) q[2];
sx q[2];
rz(-1.1897693) q[2];
sx q[2];
rz(1.5059901) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.59830211) q[1];
sx q[1];
rz(-2.1232228) q[1];
sx q[1];
rz(-1.8364947) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71535625) q[3];
sx q[3];
rz(-1.1872059) q[3];
sx q[3];
rz(0.54566783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.7168768) q[2];
sx q[2];
rz(-1.6191142) q[2];
sx q[2];
rz(2.2628722) q[2];
rz(-2.1667513) q[3];
sx q[3];
rz(-1.2474371) q[3];
sx q[3];
rz(1.6802906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.2002624) q[0];
sx q[0];
rz(-0.55568475) q[0];
sx q[0];
rz(-2.9479807) q[0];
rz(3.0337453) q[1];
sx q[1];
rz(-2.4272608) q[1];
sx q[1];
rz(-2.329619) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40859544) q[0];
sx q[0];
rz(-2.672086) q[0];
sx q[0];
rz(-1.8159189) q[0];
x q[1];
rz(2.1312563) q[2];
sx q[2];
rz(-2.3046846) q[2];
sx q[2];
rz(0.79272717) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0206611) q[1];
sx q[1];
rz(-1.1783847) q[1];
sx q[1];
rz(0.25041469) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8609686) q[3];
sx q[3];
rz(-2.6860552) q[3];
sx q[3];
rz(0.83802569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0783483) q[2];
sx q[2];
rz(-1.8599267) q[2];
sx q[2];
rz(3.1380624) q[2];
rz(-2.6053536) q[3];
sx q[3];
rz(-2.8865774) q[3];
sx q[3];
rz(-0.53853881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8435159) q[0];
sx q[0];
rz(-1.2388107) q[0];
sx q[0];
rz(2.3115944) q[0];
rz(1.4044546) q[1];
sx q[1];
rz(-2.2300215) q[1];
sx q[1];
rz(-2.8703168) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41554444) q[0];
sx q[0];
rz(-1.3726639) q[0];
sx q[0];
rz(1.8285486) q[0];
x q[1];
rz(1.6479332) q[2];
sx q[2];
rz(-2.234314) q[2];
sx q[2];
rz(-1.5653939) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8408056) q[1];
sx q[1];
rz(-1.5432165) q[1];
sx q[1];
rz(2.7843529) q[1];
rz(1.7286519) q[3];
sx q[3];
rz(-2.3612771) q[3];
sx q[3];
rz(-1.8027959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6259049) q[2];
sx q[2];
rz(-0.38334623) q[2];
sx q[2];
rz(-0.77451998) q[2];
rz(-2.1382616) q[3];
sx q[3];
rz(-2.7480875) q[3];
sx q[3];
rz(2.5956608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14458732) q[0];
sx q[0];
rz(-2.044675) q[0];
sx q[0];
rz(-1.0928094) q[0];
rz(-1.3358759) q[1];
sx q[1];
rz(-2.3922258) q[1];
sx q[1];
rz(-0.54378477) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0896801) q[0];
sx q[0];
rz(-0.26482115) q[0];
sx q[0];
rz(1.7811243) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4630578) q[2];
sx q[2];
rz(-1.4894052) q[2];
sx q[2];
rz(-0.74548364) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1256225) q[1];
sx q[1];
rz(-1.984375) q[1];
sx q[1];
rz(0.57012635) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5690992) q[3];
sx q[3];
rz(-2.5212396) q[3];
sx q[3];
rz(1.9242632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.233923) q[2];
sx q[2];
rz(-0.69207865) q[2];
sx q[2];
rz(2.862759) q[2];
rz(0.28576609) q[3];
sx q[3];
rz(-2.2279492) q[3];
sx q[3];
rz(2.7260776) q[3];
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
rz(1.435442) q[0];
sx q[0];
rz(-2.4984062) q[0];
sx q[0];
rz(-1.0253133) q[0];
rz(3.0030491) q[1];
sx q[1];
rz(-1.5899315) q[1];
sx q[1];
rz(-2.9655546) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0781949) q[0];
sx q[0];
rz(-1.0559582) q[0];
sx q[0];
rz(1.8652099) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7819809) q[2];
sx q[2];
rz(-1.7533025) q[2];
sx q[2];
rz(-0.75342853) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.61689395) q[1];
sx q[1];
rz(-2.0512676) q[1];
sx q[1];
rz(2.8887649) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.17060664) q[3];
sx q[3];
rz(-1.9904558) q[3];
sx q[3];
rz(2.4158286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.90870086) q[2];
sx q[2];
rz(-2.7140534) q[2];
sx q[2];
rz(0.78988451) q[2];
rz(0.32957736) q[3];
sx q[3];
rz(-1.7695844) q[3];
sx q[3];
rz(2.4901539) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4306507) q[0];
sx q[0];
rz(-0.21345226) q[0];
sx q[0];
rz(0.99367225) q[0];
rz(-2.1395394) q[1];
sx q[1];
rz(-2.1498945) q[1];
sx q[1];
rz(-1.4659945) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1260557) q[0];
sx q[0];
rz(-1.7483341) q[0];
sx q[0];
rz(2.4734495) q[0];
x q[1];
rz(-0.11498954) q[2];
sx q[2];
rz(-1.1144979) q[2];
sx q[2];
rz(-1.4062238) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0198145) q[1];
sx q[1];
rz(-1.3899068) q[1];
sx q[1];
rz(1.719285) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3099019) q[3];
sx q[3];
rz(-1.9134095) q[3];
sx q[3];
rz(-2.6765842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7953636) q[2];
sx q[2];
rz(-2.2987404) q[2];
sx q[2];
rz(-0.25974926) q[2];
rz(0.78884697) q[3];
sx q[3];
rz(-0.84281054) q[3];
sx q[3];
rz(1.747725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9144834) q[0];
sx q[0];
rz(-1.2863343) q[0];
sx q[0];
rz(2.1504543) q[0];
rz(1.3369075) q[1];
sx q[1];
rz(-0.82695812) q[1];
sx q[1];
rz(1.7327259) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1240272) q[0];
sx q[0];
rz(-2.6691737) q[0];
sx q[0];
rz(0.25644619) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7193871) q[2];
sx q[2];
rz(-0.5432866) q[2];
sx q[2];
rz(-3.0843184) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.91193491) q[1];
sx q[1];
rz(-1.2951944) q[1];
sx q[1];
rz(0.36987856) q[1];
rz(-1.7400319) q[3];
sx q[3];
rz(-0.90247655) q[3];
sx q[3];
rz(0.066370666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1651429) q[2];
sx q[2];
rz(-1.4850478) q[2];
sx q[2];
rz(-2.5531947) q[2];
rz(-0.31383651) q[3];
sx q[3];
rz(-2.2321759) q[3];
sx q[3];
rz(-0.61217827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9070076) q[0];
sx q[0];
rz(-2.0186277) q[0];
sx q[0];
rz(3.0871952) q[0];
rz(-0.76039487) q[1];
sx q[1];
rz(-2.1957928) q[1];
sx q[1];
rz(1.0992345) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96241073) q[0];
sx q[0];
rz(-2.1405309) q[0];
sx q[0];
rz(0.197844) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3620891) q[2];
sx q[2];
rz(-0.94346775) q[2];
sx q[2];
rz(1.9169538) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1832378) q[1];
sx q[1];
rz(-2.3545165) q[1];
sx q[1];
rz(1.3233722) q[1];
rz(-pi) q[2];
x q[2];
rz(1.681692) q[3];
sx q[3];
rz(-1.520021) q[3];
sx q[3];
rz(-2.5230356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.86193639) q[2];
sx q[2];
rz(-1.6013689) q[2];
sx q[2];
rz(-0.81735617) q[2];
rz(-3.127408) q[3];
sx q[3];
rz(-2.7458906) q[3];
sx q[3];
rz(1.0891677) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8297183) q[0];
sx q[0];
rz(-1.051396) q[0];
sx q[0];
rz(0.22148393) q[0];
rz(0.43137506) q[1];
sx q[1];
rz(-0.90310496) q[1];
sx q[1];
rz(1.0198786) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7542091) q[0];
sx q[0];
rz(-1.6176244) q[0];
sx q[0];
rz(1.55051) q[0];
x q[1];
rz(0.47616565) q[2];
sx q[2];
rz(-0.22606255) q[2];
sx q[2];
rz(3.0788596) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4012667) q[1];
sx q[1];
rz(-1.6049084) q[1];
sx q[1];
rz(-1.2383862) q[1];
rz(-0.13615578) q[3];
sx q[3];
rz(-1.7014967) q[3];
sx q[3];
rz(-3.1112373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.43883103) q[2];
sx q[2];
rz(-1.9352545) q[2];
sx q[2];
rz(-1.2552525) q[2];
rz(-3.1158279) q[3];
sx q[3];
rz(-1.3296826) q[3];
sx q[3];
rz(-0.51938957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7230566) q[0];
sx q[0];
rz(-2.3522455) q[0];
sx q[0];
rz(1.1242207) q[0];
rz(2.8070246) q[1];
sx q[1];
rz(-2.6523013) q[1];
sx q[1];
rz(2.1025067) q[1];
rz(-0.4998906) q[2];
sx q[2];
rz(-0.51639787) q[2];
sx q[2];
rz(-2.7698386) q[2];
rz(-2.6565537) q[3];
sx q[3];
rz(-2.1342604) q[3];
sx q[3];
rz(0.66935183) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
