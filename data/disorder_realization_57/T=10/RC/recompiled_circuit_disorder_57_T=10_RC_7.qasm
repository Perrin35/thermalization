OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6141619) q[0];
sx q[0];
rz(-0.80978137) q[0];
sx q[0];
rz(-0.53139395) q[0];
rz(-2.9040789) q[1];
sx q[1];
rz(-1.3637435) q[1];
sx q[1];
rz(-1.9385424) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7145568) q[0];
sx q[0];
rz(-1.3792975) q[0];
sx q[0];
rz(-1.7869851) q[0];
rz(-pi) q[1];
rz(0.25469829) q[2];
sx q[2];
rz(-1.7809976) q[2];
sx q[2];
rz(-2.3496698) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4431475) q[1];
sx q[1];
rz(-1.6277715) q[1];
sx q[1];
rz(-0.2351825) q[1];
x q[2];
rz(2.8486512) q[3];
sx q[3];
rz(-1.5353234) q[3];
sx q[3];
rz(2.5778511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.27665859) q[2];
sx q[2];
rz(-2.0099137) q[2];
sx q[2];
rz(2.9425088) q[2];
rz(1.52786) q[3];
sx q[3];
rz(-2.644643) q[3];
sx q[3];
rz(2.4075107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99877015) q[0];
sx q[0];
rz(-0.12717371) q[0];
sx q[0];
rz(0.81940991) q[0];
rz(-0.2858513) q[1];
sx q[1];
rz(-2.2276623) q[1];
sx q[1];
rz(-1.2664638) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0985384) q[0];
sx q[0];
rz(-0.96999723) q[0];
sx q[0];
rz(-1.1061125) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1663213) q[2];
sx q[2];
rz(-1.3504488) q[2];
sx q[2];
rz(-1.599556) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.61050615) q[1];
sx q[1];
rz(-0.55263457) q[1];
sx q[1];
rz(-2.3081739) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5655079) q[3];
sx q[3];
rz(-1.1844716) q[3];
sx q[3];
rz(1.9379804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1566029) q[2];
sx q[2];
rz(-1.2133602) q[2];
sx q[2];
rz(-1.161969) q[2];
rz(0.087163838) q[3];
sx q[3];
rz(-1.6379084) q[3];
sx q[3];
rz(-0.073908977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53482985) q[0];
sx q[0];
rz(-1.1207026) q[0];
sx q[0];
rz(2.8379922) q[0];
rz(-1.3820232) q[1];
sx q[1];
rz(-1.2761812) q[1];
sx q[1];
rz(0.64750013) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2979463) q[0];
sx q[0];
rz(-0.93695153) q[0];
sx q[0];
rz(-1.514939) q[0];
rz(-pi) q[1];
rz(-0.71543773) q[2];
sx q[2];
rz(-1.9999265) q[2];
sx q[2];
rz(1.5646343) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4422569) q[1];
sx q[1];
rz(-1.2067716) q[1];
sx q[1];
rz(-2.060021) q[1];
x q[2];
rz(0.89214274) q[3];
sx q[3];
rz(-2.3193079) q[3];
sx q[3];
rz(3.126614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2993762) q[2];
sx q[2];
rz(-2.3589578) q[2];
sx q[2];
rz(-1.5608609) q[2];
rz(0.98172274) q[3];
sx q[3];
rz(-1.862062) q[3];
sx q[3];
rz(0.64341199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.9075539) q[0];
sx q[0];
rz(-0.83775318) q[0];
sx q[0];
rz(-2.2696944) q[0];
rz(-0.38726989) q[1];
sx q[1];
rz(-0.64278066) q[1];
sx q[1];
rz(2.8505468) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1211348) q[0];
sx q[0];
rz(-2.1452248) q[0];
sx q[0];
rz(-2.8681884) q[0];
rz(0.83424763) q[2];
sx q[2];
rz(-0.4220037) q[2];
sx q[2];
rz(-0.01854245) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.93451553) q[1];
sx q[1];
rz(-2.0743437) q[1];
sx q[1];
rz(0.7657004) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8157186) q[3];
sx q[3];
rz(-0.59026736) q[3];
sx q[3];
rz(-3.0274689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7724093) q[2];
sx q[2];
rz(-0.79493752) q[2];
sx q[2];
rz(-0.54405653) q[2];
rz(2.6323075) q[3];
sx q[3];
rz(-1.0638758) q[3];
sx q[3];
rz(-0.86597401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0020224) q[0];
sx q[0];
rz(-1.0263356) q[0];
sx q[0];
rz(-2.496526) q[0];
rz(0.6257261) q[1];
sx q[1];
rz(-1.9179683) q[1];
sx q[1];
rz(0.63017875) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7932927) q[0];
sx q[0];
rz(-2.5508946) q[0];
sx q[0];
rz(-0.22038711) q[0];
x q[1];
rz(-2.7475287) q[2];
sx q[2];
rz(-2.2820018) q[2];
sx q[2];
rz(-1.1472536) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.1944151) q[1];
sx q[1];
rz(-1.0395323) q[1];
sx q[1];
rz(2.8814425) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.418525) q[3];
sx q[3];
rz(-2.1232743) q[3];
sx q[3];
rz(3.1395562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3877635) q[2];
sx q[2];
rz(-1.3239653) q[2];
sx q[2];
rz(1.7774263) q[2];
rz(-1.8917313) q[3];
sx q[3];
rz(-1.2156237) q[3];
sx q[3];
rz(2.2201339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.0548148) q[0];
sx q[0];
rz(-0.58371109) q[0];
sx q[0];
rz(-0.71682799) q[0];
rz(2.7038799) q[1];
sx q[1];
rz(-2.4611459) q[1];
sx q[1];
rz(-0.81370083) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2700304) q[0];
sx q[0];
rz(-0.22261482) q[0];
sx q[0];
rz(1.2827669) q[0];
x q[1];
rz(2.4196163) q[2];
sx q[2];
rz(-2.0275896) q[2];
sx q[2];
rz(-1.7320088) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3354213) q[1];
sx q[1];
rz(-0.85039447) q[1];
sx q[1];
rz(-0.45067388) q[1];
x q[2];
rz(-0.30721174) q[3];
sx q[3];
rz(-1.2324411) q[3];
sx q[3];
rz(-3.0708145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2563236) q[2];
sx q[2];
rz(-2.7041114) q[2];
sx q[2];
rz(-1.9539333) q[2];
rz(-0.93196431) q[3];
sx q[3];
rz(-2.0075802) q[3];
sx q[3];
rz(-0.37117547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1535783) q[0];
sx q[0];
rz(-2.1126641) q[0];
sx q[0];
rz(0.6860835) q[0];
rz(1.5230806) q[1];
sx q[1];
rz(-0.97421092) q[1];
sx q[1];
rz(0.66326052) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8406869) q[0];
sx q[0];
rz(-1.5947043) q[0];
sx q[0];
rz(2.2211214) q[0];
rz(-pi) q[1];
rz(-1.3528321) q[2];
sx q[2];
rz(-2.1969165) q[2];
sx q[2];
rz(-1.8916212) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.85147334) q[1];
sx q[1];
rz(-2.0514601) q[1];
sx q[1];
rz(-2.9139247) q[1];
rz(-pi) q[2];
rz(0.68526666) q[3];
sx q[3];
rz(-2.4589834) q[3];
sx q[3];
rz(0.21041378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.474581) q[2];
sx q[2];
rz(-2.3040999) q[2];
sx q[2];
rz(-0.015080301) q[2];
rz(-2.1298501) q[3];
sx q[3];
rz(-2.0254617) q[3];
sx q[3];
rz(2.570178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0886154) q[0];
sx q[0];
rz(-0.73260728) q[0];
sx q[0];
rz(0.10738871) q[0];
rz(-0.30474162) q[1];
sx q[1];
rz(-1.3184897) q[1];
sx q[1];
rz(2.6838141) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7055571) q[0];
sx q[0];
rz(-1.8920664) q[0];
sx q[0];
rz(0.9106439) q[0];
rz(-pi) q[1];
rz(0.56577487) q[2];
sx q[2];
rz(-2.3104295) q[2];
sx q[2];
rz(-2.4772252) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1220408) q[1];
sx q[1];
rz(-1.8035839) q[1];
sx q[1];
rz(1.9368534) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5910208) q[3];
sx q[3];
rz(-0.62567657) q[3];
sx q[3];
rz(-1.613137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7765939) q[2];
sx q[2];
rz(-1.685131) q[2];
sx q[2];
rz(1.7101074) q[2];
rz(-1.4252023) q[3];
sx q[3];
rz(-2.5148354) q[3];
sx q[3];
rz(-1.2861929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15437056) q[0];
sx q[0];
rz(-0.4796589) q[0];
sx q[0];
rz(-2.1881058) q[0];
rz(1.3257239) q[1];
sx q[1];
rz(-2.6486501) q[1];
sx q[1];
rz(0.58473933) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3922694) q[0];
sx q[0];
rz(-1.4587147) q[0];
sx q[0];
rz(-0.49082503) q[0];
x q[1];
rz(1.2051177) q[2];
sx q[2];
rz(-1.056991) q[2];
sx q[2];
rz(2.6363381) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4084088) q[1];
sx q[1];
rz(-1.5015263) q[1];
sx q[1];
rz(2.7015711) q[1];
x q[2];
rz(2.9986266) q[3];
sx q[3];
rz(-1.3060095) q[3];
sx q[3];
rz(2.9395482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8087625) q[2];
sx q[2];
rz(-2.3213883) q[2];
sx q[2];
rz(-0.27627036) q[2];
rz(2.590498) q[3];
sx q[3];
rz(-1.3994183) q[3];
sx q[3];
rz(0.38366693) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0693414) q[0];
sx q[0];
rz(-2.6273917) q[0];
sx q[0];
rz(0.65548354) q[0];
rz(-1.9845225) q[1];
sx q[1];
rz(-1.3815222) q[1];
sx q[1];
rz(-1.5225333) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64064202) q[0];
sx q[0];
rz(-2.2146533) q[0];
sx q[0];
rz(1.9615016) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40752878) q[2];
sx q[2];
rz(-2.4113843) q[2];
sx q[2];
rz(-0.66116316) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.22776991) q[1];
sx q[1];
rz(-1.4064944) q[1];
sx q[1];
rz(-2.5096202) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1679681) q[3];
sx q[3];
rz(-1.9625003) q[3];
sx q[3];
rz(-0.21564461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.15554252) q[2];
sx q[2];
rz(-2.7337998) q[2];
sx q[2];
rz(-0.84214169) q[2];
rz(-1.6048253) q[3];
sx q[3];
rz(-1.7611046) q[3];
sx q[3];
rz(-3.1150637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2355272) q[0];
sx q[0];
rz(-1.1885831) q[0];
sx q[0];
rz(2.6901235) q[0];
rz(-0.72262598) q[1];
sx q[1];
rz(-2.2650748) q[1];
sx q[1];
rz(1.5320019) q[1];
rz(-3.1149574) q[2];
sx q[2];
rz(-1.3616818) q[2];
sx q[2];
rz(1.5875265) q[2];
rz(1.4747254) q[3];
sx q[3];
rz(-0.69312743) q[3];
sx q[3];
rz(3.0467924) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];