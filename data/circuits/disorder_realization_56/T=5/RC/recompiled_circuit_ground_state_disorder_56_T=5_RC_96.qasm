OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.11256448) q[0];
sx q[0];
rz(-1.636314) q[0];
sx q[0];
rz(0.78328744) q[0];
rz(-2.0593491) q[1];
sx q[1];
rz(-1.6545656) q[1];
sx q[1];
rz(2.3496871) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91814268) q[0];
sx q[0];
rz(-2.8200216) q[0];
sx q[0];
rz(-2.5173223) q[0];
rz(1.804084) q[2];
sx q[2];
rz(-0.15809862) q[2];
sx q[2];
rz(2.5776517) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.72630771) q[1];
sx q[1];
rz(-1.0276762) q[1];
sx q[1];
rz(-1.607479) q[1];
rz(-pi) q[2];
rz(0.79444076) q[3];
sx q[3];
rz(-1.1367961) q[3];
sx q[3];
rz(1.3484914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.59114328) q[2];
sx q[2];
rz(-2.998896) q[2];
sx q[2];
rz(0.099451065) q[2];
rz(2.8948696) q[3];
sx q[3];
rz(-1.5244502) q[3];
sx q[3];
rz(2.7439086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0196911) q[0];
sx q[0];
rz(-0.97625232) q[0];
sx q[0];
rz(-0.9285399) q[0];
rz(-1.4506725) q[1];
sx q[1];
rz(-1.848315) q[1];
sx q[1];
rz(2.4931152) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47005338) q[0];
sx q[0];
rz(-0.82255615) q[0];
sx q[0];
rz(-0.7732735) q[0];
rz(0.5669539) q[2];
sx q[2];
rz(-1.840072) q[2];
sx q[2];
rz(-0.93967162) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8806421) q[1];
sx q[1];
rz(-1.7127258) q[1];
sx q[1];
rz(-1.5550343) q[1];
rz(1.2821372) q[3];
sx q[3];
rz(-1.4573754) q[3];
sx q[3];
rz(2.3577549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6917981) q[2];
sx q[2];
rz(-0.75642502) q[2];
sx q[2];
rz(0.85255426) q[2];
rz(2.2954588) q[3];
sx q[3];
rz(-1.6328014) q[3];
sx q[3];
rz(2.1107296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6344675) q[0];
sx q[0];
rz(-1.2256624) q[0];
sx q[0];
rz(-2.0563828) q[0];
rz(-2.197544) q[1];
sx q[1];
rz(-1.4986821) q[1];
sx q[1];
rz(-1.5012213) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2417898) q[0];
sx q[0];
rz(-2.9405624) q[0];
sx q[0];
rz(0.19636671) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3365715) q[2];
sx q[2];
rz(-1.7570436) q[2];
sx q[2];
rz(-0.062092394) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5021421) q[1];
sx q[1];
rz(-2.1413364) q[1];
sx q[1];
rz(-3.023663) q[1];
rz(-pi) q[2];
rz(1.3543868) q[3];
sx q[3];
rz(-1.7248271) q[3];
sx q[3];
rz(2.0311525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40470716) q[2];
sx q[2];
rz(-0.53386226) q[2];
sx q[2];
rz(-2.286818) q[2];
rz(0.9084304) q[3];
sx q[3];
rz(-1.5769438) q[3];
sx q[3];
rz(2.0250208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92723769) q[0];
sx q[0];
rz(-2.7404009) q[0];
sx q[0];
rz(-1.1965363) q[0];
rz(1.475097) q[1];
sx q[1];
rz(-0.42088446) q[1];
sx q[1];
rz(1.9570785) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0752206) q[0];
sx q[0];
rz(-2.3983404) q[0];
sx q[0];
rz(0.35510285) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3802211) q[2];
sx q[2];
rz(-2.3146176) q[2];
sx q[2];
rz(2.8883052) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2503916) q[1];
sx q[1];
rz(-1.7970835) q[1];
sx q[1];
rz(2.0865284) q[1];
x q[2];
rz(0.58301894) q[3];
sx q[3];
rz(-2.8172917) q[3];
sx q[3];
rz(0.7079269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.92675942) q[2];
sx q[2];
rz(-0.49761179) q[2];
sx q[2];
rz(1.8801749) q[2];
rz(2.7091806) q[3];
sx q[3];
rz(-1.132248) q[3];
sx q[3];
rz(2.6692218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.072902) q[0];
sx q[0];
rz(-1.9111159) q[0];
sx q[0];
rz(3.1121837) q[0];
rz(2.3853761) q[1];
sx q[1];
rz(-0.58719802) q[1];
sx q[1];
rz(2.3888033) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29460337) q[0];
sx q[0];
rz(-1.5651817) q[0];
sx q[0];
rz(-3.1412197) q[0];
rz(-pi) q[1];
rz(-1.2529066) q[2];
sx q[2];
rz(-0.28159062) q[2];
sx q[2];
rz(-1.3791305) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.50548762) q[1];
sx q[1];
rz(-1.6754181) q[1];
sx q[1];
rz(-0.57144798) q[1];
rz(-pi) q[2];
rz(0.54075586) q[3];
sx q[3];
rz(-1.1782559) q[3];
sx q[3];
rz(2.768571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.14043643) q[2];
sx q[2];
rz(-1.6480646) q[2];
sx q[2];
rz(2.746554) q[2];
rz(-0.47518528) q[3];
sx q[3];
rz(-1.7893712) q[3];
sx q[3];
rz(0.37676677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7659371) q[0];
sx q[0];
rz(-0.7203311) q[0];
sx q[0];
rz(-2.1790867) q[0];
rz(0.51482254) q[1];
sx q[1];
rz(-1.6074901) q[1];
sx q[1];
rz(2.8033676) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0926577) q[0];
sx q[0];
rz(-0.83507628) q[0];
sx q[0];
rz(-2.8754763) q[0];
rz(-pi) q[1];
rz(0.2372431) q[2];
sx q[2];
rz(-1.8117935) q[2];
sx q[2];
rz(-2.4324696) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.92300999) q[1];
sx q[1];
rz(-1.5708231) q[1];
sx q[1];
rz(1.567651) q[1];
rz(1.4184444) q[3];
sx q[3];
rz(-0.42644106) q[3];
sx q[3];
rz(1.8352729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6414791) q[2];
sx q[2];
rz(-1.9794455) q[2];
sx q[2];
rz(2.5872453) q[2];
rz(-2.2902299) q[3];
sx q[3];
rz(-0.35836372) q[3];
sx q[3];
rz(2.5135777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8165269) q[0];
sx q[0];
rz(-0.61282235) q[0];
sx q[0];
rz(-2.391173) q[0];
rz(0.57506192) q[1];
sx q[1];
rz(-1.776418) q[1];
sx q[1];
rz(-2.4837928) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0477792) q[0];
sx q[0];
rz(-1.1623315) q[0];
sx q[0];
rz(-0.7042709) q[0];
rz(1.6148241) q[2];
sx q[2];
rz(-0.86404534) q[2];
sx q[2];
rz(2.6296774) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.21322528) q[1];
sx q[1];
rz(-1.2847273) q[1];
sx q[1];
rz(3.1242773) q[1];
rz(-pi) q[2];
rz(-1.1893473) q[3];
sx q[3];
rz(-1.369841) q[3];
sx q[3];
rz(0.23924669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.647992) q[2];
sx q[2];
rz(-1.0085663) q[2];
sx q[2];
rz(1.4823401) q[2];
rz(0.12065398) q[3];
sx q[3];
rz(-1.3841265) q[3];
sx q[3];
rz(2.306126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3487314) q[0];
sx q[0];
rz(-0.86935765) q[0];
sx q[0];
rz(-2.8435775) q[0];
rz(-2.1030078) q[1];
sx q[1];
rz(-2.5972001) q[1];
sx q[1];
rz(0.15377741) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6385495) q[0];
sx q[0];
rz(-2.8766603) q[0];
sx q[0];
rz(2.8517836) q[0];
x q[1];
rz(2.2258899) q[2];
sx q[2];
rz(-2.2580552) q[2];
sx q[2];
rz(-2.5926431) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7419305) q[1];
sx q[1];
rz(-2.7832418) q[1];
sx q[1];
rz(-2.8835924) q[1];
rz(-pi) q[2];
x q[2];
rz(0.12184398) q[3];
sx q[3];
rz(-1.9949556) q[3];
sx q[3];
rz(0.76790732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.095857233) q[2];
sx q[2];
rz(-0.46892527) q[2];
sx q[2];
rz(-1.9528961) q[2];
rz(-1.3304322) q[3];
sx q[3];
rz(-1.1878139) q[3];
sx q[3];
rz(-2.4082898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7981912) q[0];
sx q[0];
rz(-0.60438406) q[0];
sx q[0];
rz(-1.4991624) q[0];
rz(-2.5921953) q[1];
sx q[1];
rz(-1.5769985) q[1];
sx q[1];
rz(-0.25064358) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0325913) q[0];
sx q[0];
rz(-1.8792986) q[0];
sx q[0];
rz(0.89621131) q[0];
rz(-pi) q[1];
x q[1];
rz(0.05392404) q[2];
sx q[2];
rz(-2.8802875) q[2];
sx q[2];
rz(-0.037029412) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0359866) q[1];
sx q[1];
rz(-0.79933724) q[1];
sx q[1];
rz(1.4738333) q[1];
x q[2];
rz(0.57119675) q[3];
sx q[3];
rz(-0.2303752) q[3];
sx q[3];
rz(-0.26709589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9833019) q[2];
sx q[2];
rz(-3.0088708) q[2];
sx q[2];
rz(-1.0164725) q[2];
rz(0.086183444) q[3];
sx q[3];
rz(-2.1250171) q[3];
sx q[3];
rz(-2.3763954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30855274) q[0];
sx q[0];
rz(-0.85423952) q[0];
sx q[0];
rz(2.7225851) q[0];
rz(2.6047756) q[1];
sx q[1];
rz(-0.62428004) q[1];
sx q[1];
rz(0.73582617) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9402855) q[0];
sx q[0];
rz(-0.92865151) q[0];
sx q[0];
rz(-1.3298195) q[0];
x q[1];
rz(2.4707253) q[2];
sx q[2];
rz(-2.3773674) q[2];
sx q[2];
rz(-2.4799181) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.23259737) q[1];
sx q[1];
rz(-0.67282644) q[1];
sx q[1];
rz(-2.7907967) q[1];
rz(-pi) q[2];
rz(-2.5352468) q[3];
sx q[3];
rz(-2.5744573) q[3];
sx q[3];
rz(0.66886574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.95911038) q[2];
sx q[2];
rz(-0.82254326) q[2];
sx q[2];
rz(-0.62270069) q[2];
rz(-0.73838082) q[3];
sx q[3];
rz(-1.9564956) q[3];
sx q[3];
rz(0.27899376) q[3];
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
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5644792) q[0];
sx q[0];
rz(-2.8688685) q[0];
sx q[0];
rz(-2.175749) q[0];
rz(-0.64361698) q[1];
sx q[1];
rz(-1.2871965) q[1];
sx q[1];
rz(1.0866477) q[1];
rz(1.9401445) q[2];
sx q[2];
rz(-1.7864173) q[2];
sx q[2];
rz(-1.3997072) q[2];
rz(0.24091992) q[3];
sx q[3];
rz(-1.708247) q[3];
sx q[3];
rz(0.74549992) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
