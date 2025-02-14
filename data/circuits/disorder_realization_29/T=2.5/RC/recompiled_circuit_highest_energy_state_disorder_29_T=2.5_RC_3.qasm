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
rz(1.2443378) q[0];
sx q[0];
rz(-0.79255784) q[0];
sx q[0];
rz(2.8886524) q[0];
rz(-0.77415544) q[1];
sx q[1];
rz(-0.3781265) q[1];
sx q[1];
rz(0.3869431) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16432504) q[0];
sx q[0];
rz(-1.2480134) q[0];
sx q[0];
rz(2.8708145) q[0];
x q[1];
rz(-1.4841485) q[2];
sx q[2];
rz(-0.99471015) q[2];
sx q[2];
rz(-0.91743166) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.746443) q[1];
sx q[1];
rz(-0.11583318) q[1];
sx q[1];
rz(-2.3833116) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4693842) q[3];
sx q[3];
rz(-0.70939964) q[3];
sx q[3];
rz(-1.0649452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0953377) q[2];
sx q[2];
rz(-0.75901186) q[2];
sx q[2];
rz(-1.7389899) q[2];
rz(-1.4143573) q[3];
sx q[3];
rz(-0.71310133) q[3];
sx q[3];
rz(0.49899092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6801179) q[0];
sx q[0];
rz(-2.6533227) q[0];
sx q[0];
rz(-2.4312191) q[0];
rz(1.0164227) q[1];
sx q[1];
rz(-1.2998394) q[1];
sx q[1];
rz(-2.3764835) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23952627) q[0];
sx q[0];
rz(-1.2529116) q[0];
sx q[0];
rz(-0.56196281) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76025195) q[2];
sx q[2];
rz(-2.2905802) q[2];
sx q[2];
rz(-0.30530294) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2110846) q[1];
sx q[1];
rz(-0.74568664) q[1];
sx q[1];
rz(1.0812283) q[1];
rz(-2.9917172) q[3];
sx q[3];
rz(-1.2399763) q[3];
sx q[3];
rz(-0.56043032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1171099) q[2];
sx q[2];
rz(-2.4958002) q[2];
sx q[2];
rz(0.70466858) q[2];
rz(2.9674528) q[3];
sx q[3];
rz(-0.87248674) q[3];
sx q[3];
rz(2.7599938) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0524549) q[0];
sx q[0];
rz(-1.8280886) q[0];
sx q[0];
rz(-2.254159) q[0];
rz(1.2499836) q[1];
sx q[1];
rz(-1.9307815) q[1];
sx q[1];
rz(-1.3608305) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0075390752) q[0];
sx q[0];
rz(-2.0749395) q[0];
sx q[0];
rz(-1.1372805) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7086231) q[2];
sx q[2];
rz(-1.6500452) q[2];
sx q[2];
rz(-2.2887231) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2584784) q[1];
sx q[1];
rz(-0.7683903) q[1];
sx q[1];
rz(-1.7151296) q[1];
rz(-pi) q[2];
rz(3.0988337) q[3];
sx q[3];
rz(-1.1819289) q[3];
sx q[3];
rz(-2.1270909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.36797324) q[2];
sx q[2];
rz(-2.7757288) q[2];
sx q[2];
rz(1.6486637) q[2];
rz(2.6922373) q[3];
sx q[3];
rz(-1.8466693) q[3];
sx q[3];
rz(1.2202643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.087695) q[0];
sx q[0];
rz(-0.59006417) q[0];
sx q[0];
rz(1.4087403) q[0];
rz(-0.98211163) q[1];
sx q[1];
rz(-1.5928007) q[1];
sx q[1];
rz(1.7353479) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81534751) q[0];
sx q[0];
rz(-2.8994588) q[0];
sx q[0];
rz(1.0905488) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6236308) q[2];
sx q[2];
rz(-1.5312734) q[2];
sx q[2];
rz(1.7588577) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6342369) q[1];
sx q[1];
rz(-2.0869531) q[1];
sx q[1];
rz(-2.8576351) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9191756) q[3];
sx q[3];
rz(-2.4628785) q[3];
sx q[3];
rz(1.410759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9881607) q[2];
sx q[2];
rz(-2.1250696) q[2];
sx q[2];
rz(-0.15596685) q[2];
rz(2.4912452) q[3];
sx q[3];
rz(-1.1740843) q[3];
sx q[3];
rz(3.0273738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0680189) q[0];
sx q[0];
rz(-1.8719712) q[0];
sx q[0];
rz(0.20075783) q[0];
rz(1.1772032) q[1];
sx q[1];
rz(-2.2889844) q[1];
sx q[1];
rz(-1.9815365) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2748799) q[0];
sx q[0];
rz(-1.794941) q[0];
sx q[0];
rz(0.17205162) q[0];
rz(-pi) q[1];
rz(0.11518466) q[2];
sx q[2];
rz(-1.6403926) q[2];
sx q[2];
rz(1.4032961) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0811942) q[1];
sx q[1];
rz(-1.8738998) q[1];
sx q[1];
rz(-0.36784192) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0409768) q[3];
sx q[3];
rz(-1.3966832) q[3];
sx q[3];
rz(-1.8903738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3843627) q[2];
sx q[2];
rz(-0.94567662) q[2];
sx q[2];
rz(-0.45822701) q[2];
rz(0.630817) q[3];
sx q[3];
rz(-1.9061371) q[3];
sx q[3];
rz(0.49921504) q[3];
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
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45796564) q[0];
sx q[0];
rz(-0.85245913) q[0];
sx q[0];
rz(-0.0068579554) q[0];
rz(-0.231617) q[1];
sx q[1];
rz(-1.7624785) q[1];
sx q[1];
rz(-2.0793656) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7201014) q[0];
sx q[0];
rz(-1.4473697) q[0];
sx q[0];
rz(-0.31436679) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3590603) q[2];
sx q[2];
rz(-2.0604302) q[2];
sx q[2];
rz(-0.28070606) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0343218) q[1];
sx q[1];
rz(-2.2244029) q[1];
sx q[1];
rz(-1.2493709) q[1];
x q[2];
rz(-0.16522629) q[3];
sx q[3];
rz(-2.2731785) q[3];
sx q[3];
rz(0.14152292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.14083938) q[2];
sx q[2];
rz(-1.2522298) q[2];
sx q[2];
rz(-1.2678649) q[2];
rz(0.86152348) q[3];
sx q[3];
rz(-1.0637161) q[3];
sx q[3];
rz(-1.4388194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2086901) q[0];
sx q[0];
rz(-0.33354315) q[0];
sx q[0];
rz(-0.21555899) q[0];
rz(-2.6453099) q[1];
sx q[1];
rz(-1.1880621) q[1];
sx q[1];
rz(2.9821679) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2886141) q[0];
sx q[0];
rz(-2.5067177) q[0];
sx q[0];
rz(1.575241) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47093289) q[2];
sx q[2];
rz(-1.3976946) q[2];
sx q[2];
rz(1.7488232) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.39495271) q[1];
sx q[1];
rz(-1.3356908) q[1];
sx q[1];
rz(3.0993153) q[1];
x q[2];
rz(-2.6786118) q[3];
sx q[3];
rz(-1.0701051) q[3];
sx q[3];
rz(2.7241796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7364007) q[2];
sx q[2];
rz(-1.1242194) q[2];
sx q[2];
rz(2.6250725) q[2];
rz(1.2107595) q[3];
sx q[3];
rz(-1.4529994) q[3];
sx q[3];
rz(1.3794544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6571758) q[0];
sx q[0];
rz(-0.076198904) q[0];
sx q[0];
rz(-2.6013689) q[0];
rz(-1.6172488) q[1];
sx q[1];
rz(-2.1397619) q[1];
sx q[1];
rz(-1.2670955) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57813533) q[0];
sx q[0];
rz(-1.7892013) q[0];
sx q[0];
rz(-2.6294623) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2742496) q[2];
sx q[2];
rz(-0.92716414) q[2];
sx q[2];
rz(0.95736733) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.051991302) q[1];
sx q[1];
rz(-0.39685186) q[1];
sx q[1];
rz(-2.7605935) q[1];
x q[2];
rz(-2.0319875) q[3];
sx q[3];
rz(-2.4465585) q[3];
sx q[3];
rz(-2.024533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2544864) q[2];
sx q[2];
rz(-1.1595414) q[2];
sx q[2];
rz(-2.9642588) q[2];
rz(1.0698498) q[3];
sx q[3];
rz(-2.0900574) q[3];
sx q[3];
rz(-0.18226084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7808481) q[0];
sx q[0];
rz(-1.9875263) q[0];
sx q[0];
rz(-0.10072197) q[0];
rz(-1.992647) q[1];
sx q[1];
rz(-1.9457685) q[1];
sx q[1];
rz(-1.9675868) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35648221) q[0];
sx q[0];
rz(-2.504341) q[0];
sx q[0];
rz(-2.8204945) q[0];
x q[1];
rz(0.76624932) q[2];
sx q[2];
rz(-1.5011929) q[2];
sx q[2];
rz(-1.0468259) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8496597) q[1];
sx q[1];
rz(-1.5330452) q[1];
sx q[1];
rz(0.019456073) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10281201) q[3];
sx q[3];
rz(-2.2592756) q[3];
sx q[3];
rz(-1.4172872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9775057) q[2];
sx q[2];
rz(-1.6567433) q[2];
sx q[2];
rz(2.738415) q[2];
rz(-2.4781503) q[3];
sx q[3];
rz(-2.5622029) q[3];
sx q[3];
rz(-3.1373533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73096257) q[0];
sx q[0];
rz(-1.0799438) q[0];
sx q[0];
rz(-0.078911111) q[0];
rz(2.2321189) q[1];
sx q[1];
rz(-2.0957004) q[1];
sx q[1];
rz(-2.7452819) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8709804) q[0];
sx q[0];
rz(-1.5798693) q[0];
sx q[0];
rz(-1.5551644) q[0];
rz(1.5635919) q[2];
sx q[2];
rz(-2.3644937) q[2];
sx q[2];
rz(-0.33852984) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.24416284) q[1];
sx q[1];
rz(-1.8051935) q[1];
sx q[1];
rz(-2.6525556) q[1];
rz(-2.590606) q[3];
sx q[3];
rz(-2.1175623) q[3];
sx q[3];
rz(1.8282229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4497455) q[2];
sx q[2];
rz(-0.88374603) q[2];
sx q[2];
rz(-0.79908243) q[2];
rz(3.0633022) q[3];
sx q[3];
rz(-1.2543863) q[3];
sx q[3];
rz(-2.4962795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4770724) q[0];
sx q[0];
rz(-1.7417396) q[0];
sx q[0];
rz(-2.59792) q[0];
rz(1.1935344) q[1];
sx q[1];
rz(-1.8652893) q[1];
sx q[1];
rz(0.51020772) q[1];
rz(2.0281744) q[2];
sx q[2];
rz(-1.4174247) q[2];
sx q[2];
rz(-2.4658801) q[2];
rz(-2.4166475) q[3];
sx q[3];
rz(-1.1851539) q[3];
sx q[3];
rz(-0.44222587) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
