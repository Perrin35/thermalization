OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2064535) q[0];
sx q[0];
rz(-0.78092617) q[0];
sx q[0];
rz(-0.20679064) q[0];
rz(-2.7517154) q[1];
sx q[1];
rz(-2.0808527) q[1];
sx q[1];
rz(-2.9614255) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5500096) q[0];
sx q[0];
rz(-1.2342493) q[0];
sx q[0];
rz(-0.49206375) q[0];
x q[1];
rz(1.4397058) q[2];
sx q[2];
rz(-2.0799473) q[2];
sx q[2];
rz(-2.5703562) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.78046103) q[1];
sx q[1];
rz(-1.6447004) q[1];
sx q[1];
rz(2.2547045) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7546685) q[3];
sx q[3];
rz(-2.2443218) q[3];
sx q[3];
rz(2.9415188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7303598) q[2];
sx q[2];
rz(-0.87324548) q[2];
sx q[2];
rz(-1.2940548) q[2];
rz(-2.7358352) q[3];
sx q[3];
rz(-1.5016705) q[3];
sx q[3];
rz(0.40677235) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2186573) q[0];
sx q[0];
rz(-1.0359534) q[0];
sx q[0];
rz(-3.0157715) q[0];
rz(-2.3361092) q[1];
sx q[1];
rz(-2.3352354) q[1];
sx q[1];
rz(-1.7696101) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94700891) q[0];
sx q[0];
rz(-2.2451375) q[0];
sx q[0];
rz(-1.4916923) q[0];
rz(-pi) q[1];
x q[1];
rz(0.20184529) q[2];
sx q[2];
rz(-2.5864374) q[2];
sx q[2];
rz(1.3943878) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.068472915) q[1];
sx q[1];
rz(-1.8205376) q[1];
sx q[1];
rz(1.5239034) q[1];
rz(-pi) q[2];
rz(-2.6729229) q[3];
sx q[3];
rz(-1.7892924) q[3];
sx q[3];
rz(-0.0056497638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.25386086) q[2];
sx q[2];
rz(-2.4527206) q[2];
sx q[2];
rz(2.1453693) q[2];
rz(1.0960724) q[3];
sx q[3];
rz(-1.4527861) q[3];
sx q[3];
rz(-2.2912099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6753733) q[0];
sx q[0];
rz(-2.1815364) q[0];
sx q[0];
rz(-2.0825785) q[0];
rz(-1.9937218) q[1];
sx q[1];
rz(-0.73906001) q[1];
sx q[1];
rz(2.0770729) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7721467) q[0];
sx q[0];
rz(-1.9659974) q[0];
sx q[0];
rz(2.4482083) q[0];
x q[1];
rz(2.9840165) q[2];
sx q[2];
rz(-1.1925979) q[2];
sx q[2];
rz(-2.6842897) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3435622) q[1];
sx q[1];
rz(-2.8312632) q[1];
sx q[1];
rz(-1.6088435) q[1];
rz(-pi) q[2];
rz(1.3833984) q[3];
sx q[3];
rz(-2.4083532) q[3];
sx q[3];
rz(-0.76776615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.069783) q[2];
sx q[2];
rz(-1.9548543) q[2];
sx q[2];
rz(-0.26322571) q[2];
rz(2.0227382) q[3];
sx q[3];
rz(-1.2759195) q[3];
sx q[3];
rz(-2.3222205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6333106) q[0];
sx q[0];
rz(-3.0968956) q[0];
sx q[0];
rz(1.7472965) q[0];
rz(1.0385723) q[1];
sx q[1];
rz(-1.4515406) q[1];
sx q[1];
rz(2.0746453) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2571714) q[0];
sx q[0];
rz(-2.1002249) q[0];
sx q[0];
rz(0.04495312) q[0];
rz(-pi) q[1];
x q[1];
rz(0.28654751) q[2];
sx q[2];
rz(-0.28038014) q[2];
sx q[2];
rz(0.36975158) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2414788) q[1];
sx q[1];
rz(-0.80726868) q[1];
sx q[1];
rz(0.57358731) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0889441) q[3];
sx q[3];
rz(-1.7715766) q[3];
sx q[3];
rz(2.5364385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2944494) q[2];
sx q[2];
rz(-1.3501945) q[2];
sx q[2];
rz(2.8520544) q[2];
rz(-2.6211522) q[3];
sx q[3];
rz(-1.8013022) q[3];
sx q[3];
rz(-2.382544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6510058) q[0];
sx q[0];
rz(-0.0087954272) q[0];
sx q[0];
rz(-2.3068413) q[0];
rz(-1.2443776) q[1];
sx q[1];
rz(-1.8908187) q[1];
sx q[1];
rz(1.429819) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.017905047) q[0];
sx q[0];
rz(-0.062739685) q[0];
sx q[0];
rz(-0.040266589) q[0];
x q[1];
rz(-1.8102112) q[2];
sx q[2];
rz(-1.8481701) q[2];
sx q[2];
rz(-2.036236) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1289542) q[1];
sx q[1];
rz(-1.5476989) q[1];
sx q[1];
rz(0.79146339) q[1];
rz(0.91246446) q[3];
sx q[3];
rz(-1.4294799) q[3];
sx q[3];
rz(0.055012881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4804068) q[2];
sx q[2];
rz(-1.0057665) q[2];
sx q[2];
rz(1.2832114) q[2];
rz(-0.13218203) q[3];
sx q[3];
rz(-0.32871267) q[3];
sx q[3];
rz(0.33974084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64602393) q[0];
sx q[0];
rz(-3.0806354) q[0];
sx q[0];
rz(-2.6640889) q[0];
rz(1.6409138) q[1];
sx q[1];
rz(-1.6093107) q[1];
sx q[1];
rz(-0.57055155) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43317023) q[0];
sx q[0];
rz(-0.47009429) q[0];
sx q[0];
rz(2.2545459) q[0];
rz(-pi) q[1];
rz(-0.45285593) q[2];
sx q[2];
rz(-1.9100683) q[2];
sx q[2];
rz(-2.0138182) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.79710302) q[1];
sx q[1];
rz(-0.87901607) q[1];
sx q[1];
rz(0.89747353) q[1];
x q[2];
rz(-2.493294) q[3];
sx q[3];
rz(-1.5530506) q[3];
sx q[3];
rz(-1.363021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.70696124) q[2];
sx q[2];
rz(-2.091445) q[2];
sx q[2];
rz(2.3449507) q[2];
rz(0.32026511) q[3];
sx q[3];
rz(-1.0864778) q[3];
sx q[3];
rz(-1.8626574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0657848) q[0];
sx q[0];
rz(-1.0362352) q[0];
sx q[0];
rz(0.35807034) q[0];
rz(-0.2886731) q[1];
sx q[1];
rz(-2.6117548) q[1];
sx q[1];
rz(1.8310865) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3863556) q[0];
sx q[0];
rz(-2.1023395) q[0];
sx q[0];
rz(-1.7887572) q[0];
x q[1];
rz(-0.99408044) q[2];
sx q[2];
rz(-0.95816441) q[2];
sx q[2];
rz(-0.47555579) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3086991) q[1];
sx q[1];
rz(-0.70736865) q[1];
sx q[1];
rz(-2.4127712) q[1];
rz(0.88498022) q[3];
sx q[3];
rz(-1.0776057) q[3];
sx q[3];
rz(0.54515884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8075809) q[2];
sx q[2];
rz(-0.58471218) q[2];
sx q[2];
rz(2.1441148) q[2];
rz(2.5618662) q[3];
sx q[3];
rz(-0.72967356) q[3];
sx q[3];
rz(1.4355481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(0.26577935) q[0];
sx q[0];
rz(-1.6515323) q[0];
sx q[0];
rz(0.34061256) q[0];
rz(-1.9050725) q[1];
sx q[1];
rz(-2.3842936) q[1];
sx q[1];
rz(1.9514726) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6105881) q[0];
sx q[0];
rz(-1.6325145) q[0];
sx q[0];
rz(1.8025663) q[0];
x q[1];
rz(3.0731593) q[2];
sx q[2];
rz(-1.4901731) q[2];
sx q[2];
rz(1.8552519) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0654046) q[1];
sx q[1];
rz(-1.3082062) q[1];
sx q[1];
rz(0.29448387) q[1];
rz(-pi) q[2];
rz(2.6885919) q[3];
sx q[3];
rz(-2.0815597) q[3];
sx q[3];
rz(-0.1764899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3999346) q[2];
sx q[2];
rz(-2.2592893) q[2];
sx q[2];
rz(0.41440543) q[2];
rz(1.7587781) q[3];
sx q[3];
rz(-1.4004935) q[3];
sx q[3];
rz(-2.8167021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25093108) q[0];
sx q[0];
rz(-1.0216167) q[0];
sx q[0];
rz(-0.1517621) q[0];
rz(1.4258619) q[1];
sx q[1];
rz(-1.7646004) q[1];
sx q[1];
rz(-2.192416) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.855809) q[0];
sx q[0];
rz(-1.7109509) q[0];
sx q[0];
rz(0.44482081) q[0];
x q[1];
rz(-1.7911712) q[2];
sx q[2];
rz(-0.95249635) q[2];
sx q[2];
rz(-0.097620336) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.28724972) q[1];
sx q[1];
rz(-2.2209475) q[1];
sx q[1];
rz(1.8734422) q[1];
rz(-2.6563431) q[3];
sx q[3];
rz(-1.0356379) q[3];
sx q[3];
rz(-0.29569611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.40733797) q[2];
sx q[2];
rz(-2.7463425) q[2];
sx q[2];
rz(-0.43241832) q[2];
rz(1.5405103) q[3];
sx q[3];
rz(-1.5099022) q[3];
sx q[3];
rz(-0.66175118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0712873) q[0];
sx q[0];
rz(-2.8685331) q[0];
sx q[0];
rz(-0.37316698) q[0];
rz(2.7846653) q[1];
sx q[1];
rz(-2.280805) q[1];
sx q[1];
rz(-0.7235136) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.824677) q[0];
sx q[0];
rz(-2.4002889) q[0];
sx q[0];
rz(2.9497428) q[0];
rz(1.4805484) q[2];
sx q[2];
rz(-0.89333488) q[2];
sx q[2];
rz(1.8842763) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3011303) q[1];
sx q[1];
rz(-2.3311619) q[1];
sx q[1];
rz(-2.7644972) q[1];
x q[2];
rz(-1.5756597) q[3];
sx q[3];
rz(-1.8412207) q[3];
sx q[3];
rz(1.5750386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2202806) q[2];
sx q[2];
rz(-1.3833412) q[2];
sx q[2];
rz(-2.5881361) q[2];
rz(0.71183318) q[3];
sx q[3];
rz(-0.40829855) q[3];
sx q[3];
rz(-1.0420943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2198467) q[0];
sx q[0];
rz(-2.5384359) q[0];
sx q[0];
rz(0.13721101) q[0];
rz(2.8593821) q[1];
sx q[1];
rz(-2.3574775) q[1];
sx q[1];
rz(-2.2025253) q[1];
rz(2.0201335) q[2];
sx q[2];
rz(-1.1171787) q[2];
sx q[2];
rz(1.6906307) q[2];
rz(-1.8525193) q[3];
sx q[3];
rz(-1.9586133) q[3];
sx q[3];
rz(-1.9546399) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
