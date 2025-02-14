OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.29326987) q[0];
sx q[0];
rz(3.3942437) q[0];
sx q[0];
rz(10.630339) q[0];
rz(2.0134917) q[1];
sx q[1];
rz(-1.3265346) q[1];
sx q[1];
rz(0.74572745) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099961258) q[0];
sx q[0];
rz(-1.7929165) q[0];
sx q[0];
rz(0.061454031) q[0];
rz(-pi) q[1];
rz(-0.029736515) q[2];
sx q[2];
rz(-0.22448891) q[2];
sx q[2];
rz(2.9475074) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2499143) q[1];
sx q[1];
rz(-0.39977705) q[1];
sx q[1];
rz(-1.0067902) q[1];
rz(-pi) q[2];
rz(-0.20799237) q[3];
sx q[3];
rz(-2.0361414) q[3];
sx q[3];
rz(-0.057010827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.91743177) q[2];
sx q[2];
rz(-1.5956343) q[2];
sx q[2];
rz(1.6708299) q[2];
rz(-1.6338232) q[3];
sx q[3];
rz(-3.1247415) q[3];
sx q[3];
rz(-2.2182218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5725752) q[0];
sx q[0];
rz(-1.1976396) q[0];
sx q[0];
rz(-1.5684599) q[0];
rz(-2.9720427) q[1];
sx q[1];
rz(-0.1145656) q[1];
sx q[1];
rz(0.13470185) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7260774) q[0];
sx q[0];
rz(-1.2252136) q[0];
sx q[0];
rz(0.4536566) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7472166) q[2];
sx q[2];
rz(-3.0718832) q[2];
sx q[2];
rz(1.7084054) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.33011064) q[1];
sx q[1];
rz(-0.54381424) q[1];
sx q[1];
rz(0.61119975) q[1];
rz(1.7487583) q[3];
sx q[3];
rz(-1.6497532) q[3];
sx q[3];
rz(0.7782026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1255101) q[2];
sx q[2];
rz(-1.4849911) q[2];
sx q[2];
rz(-0.15277319) q[2];
rz(1.7759391) q[3];
sx q[3];
rz(-0.036866166) q[3];
sx q[3];
rz(2.9735145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.633054) q[0];
sx q[0];
rz(-0.80798739) q[0];
sx q[0];
rz(-0.48164865) q[0];
rz(-2.956849) q[1];
sx q[1];
rz(-1.7748723) q[1];
sx q[1];
rz(2.1462323) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5115968) q[0];
sx q[0];
rz(-0.42034402) q[0];
sx q[0];
rz(2.1612801) q[0];
rz(-pi) q[1];
rz(-3.0929933) q[2];
sx q[2];
rz(-1.6311797) q[2];
sx q[2];
rz(0.27602613) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9650594) q[1];
sx q[1];
rz(-2.4048724) q[1];
sx q[1];
rz(0.45335575) q[1];
rz(-pi) q[2];
rz(-1.5257201) q[3];
sx q[3];
rz(-2.5451676) q[3];
sx q[3];
rz(-0.10208043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2153726) q[2];
sx q[2];
rz(-0.054457713) q[2];
sx q[2];
rz(0.064662956) q[2];
rz(1.0489382) q[3];
sx q[3];
rz(-0.026853042) q[3];
sx q[3];
rz(-1.2803199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8365086) q[0];
sx q[0];
rz(-3.0236112) q[0];
sx q[0];
rz(-0.82103658) q[0];
rz(0.07846421) q[1];
sx q[1];
rz(-1.4469701) q[1];
sx q[1];
rz(-0.99748126) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2729014) q[0];
sx q[0];
rz(-1.8857546) q[0];
sx q[0];
rz(-0.64906831) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5115963) q[2];
sx q[2];
rz(-1.6069876) q[2];
sx q[2];
rz(-0.22360392) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.32179896) q[1];
sx q[1];
rz(-2.5410278) q[1];
sx q[1];
rz(2.0482045) q[1];
rz(-pi) q[2];
rz(-1.5134345) q[3];
sx q[3];
rz(-2.0485176) q[3];
sx q[3];
rz(-0.78335947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0923826) q[2];
sx q[2];
rz(-3.1123078) q[2];
sx q[2];
rz(-1.9878261) q[2];
rz(-0.18618259) q[3];
sx q[3];
rz(-3.0598873) q[3];
sx q[3];
rz(-0.33128273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1368197) q[0];
sx q[0];
rz(-0.81757075) q[0];
sx q[0];
rz(1.9983043) q[0];
rz(1.1413057) q[1];
sx q[1];
rz(-2.3316796) q[1];
sx q[1];
rz(2.6236261) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85803131) q[0];
sx q[0];
rz(-1.0378583) q[0];
sx q[0];
rz(0.74646797) q[0];
rz(-pi) q[1];
rz(0.00092351726) q[2];
sx q[2];
rz(-1.5858272) q[2];
sx q[2];
rz(1.3236486) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3822381) q[1];
sx q[1];
rz(-1.2975946) q[1];
sx q[1];
rz(-1.2279991) q[1];
rz(-0.2980026) q[3];
sx q[3];
rz(-1.4299222) q[3];
sx q[3];
rz(-0.44671392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3190069) q[2];
sx q[2];
rz(-2.1876882) q[2];
sx q[2];
rz(-2.8138568) q[2];
rz(1.1945126) q[3];
sx q[3];
rz(-3.0142398) q[3];
sx q[3];
rz(0.85668844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0026534) q[0];
sx q[0];
rz(-0.26873538) q[0];
sx q[0];
rz(-0.56185454) q[0];
rz(-1.6506763) q[1];
sx q[1];
rz(-1.5166538) q[1];
sx q[1];
rz(-3.0432826) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1931629) q[0];
sx q[0];
rz(-1.6062482) q[0];
sx q[0];
rz(-1.9928785) q[0];
rz(-pi) q[1];
rz(-0.19616429) q[2];
sx q[2];
rz(-3.1409396) q[2];
sx q[2];
rz(-0.12015039) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1001818) q[1];
sx q[1];
rz(-1.5242002) q[1];
sx q[1];
rz(0.46943922) q[1];
rz(1.1071854) q[3];
sx q[3];
rz(-2.0258198) q[3];
sx q[3];
rz(-0.95038271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.085122846) q[2];
sx q[2];
rz(-0.19530185) q[2];
sx q[2];
rz(-1.0657715) q[2];
rz(2.7417475) q[3];
sx q[3];
rz(-0.53374922) q[3];
sx q[3];
rz(-1.8736418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14168508) q[0];
sx q[0];
rz(-3.0596924) q[0];
sx q[0];
rz(1.4578693) q[0];
rz(-1.1204002) q[1];
sx q[1];
rz(-3.0034062) q[1];
sx q[1];
rz(-0.33946005) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0507815) q[0];
sx q[0];
rz(-3.0635186) q[0];
sx q[0];
rz(-2.9705621) q[0];
rz(-pi) q[1];
rz(2.5869114) q[2];
sx q[2];
rz(-1.5590057) q[2];
sx q[2];
rz(1.5654711) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.554018) q[1];
sx q[1];
rz(-1.6563935) q[1];
sx q[1];
rz(0.0016335131) q[1];
x q[2];
rz(-1.4711597) q[3];
sx q[3];
rz(-1.2335586) q[3];
sx q[3];
rz(-1.8189614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.3642984) q[2];
sx q[2];
rz(-0.0019145049) q[2];
sx q[2];
rz(0.36336362) q[2];
rz(1.0904788) q[3];
sx q[3];
rz(-0.57991475) q[3];
sx q[3];
rz(2.0232078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5217487) q[0];
sx q[0];
rz(-2.813297) q[0];
sx q[0];
rz(-2.2186665) q[0];
rz(-1.6587616) q[1];
sx q[1];
rz(-2.5117579) q[1];
sx q[1];
rz(0.0042075687) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9516966) q[0];
sx q[0];
rz(-1.0605264) q[0];
sx q[0];
rz(-2.4416591) q[0];
rz(-pi) q[1];
rz(1.5571874) q[2];
sx q[2];
rz(-2.7306692) q[2];
sx q[2];
rz(0.027337242) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.15531047) q[1];
sx q[1];
rz(-2.661507) q[1];
sx q[1];
rz(-1.7106777) q[1];
rz(3.0757853) q[3];
sx q[3];
rz(-2.5281457) q[3];
sx q[3];
rz(0.25019161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5795472) q[2];
sx q[2];
rz(-1.5374708) q[2];
sx q[2];
rz(-1.2008249) q[2];
rz(1.3935401) q[3];
sx q[3];
rz(-3.1378919) q[3];
sx q[3];
rz(-2.4414731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30018184) q[0];
sx q[0];
rz(-2.6926079) q[0];
sx q[0];
rz(1.2196983) q[0];
rz(1.8118743) q[1];
sx q[1];
rz(-1.1390353) q[1];
sx q[1];
rz(-2.9776998) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3541479) q[0];
sx q[0];
rz(-1.3992926) q[0];
sx q[0];
rz(-1.9308912) q[0];
rz(0.40028769) q[2];
sx q[2];
rz(-2.5185555) q[2];
sx q[2];
rz(-1.7965339) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6170609) q[1];
sx q[1];
rz(-1.6304029) q[1];
sx q[1];
rz(1.4527133) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2396062) q[3];
sx q[3];
rz(-2.1928582) q[3];
sx q[3];
rz(-0.40611503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.71358744) q[2];
sx q[2];
rz(-3.0516477) q[2];
sx q[2];
rz(-0.62291992) q[2];
rz(0.04341393) q[3];
sx q[3];
rz(-0.89689887) q[3];
sx q[3];
rz(-2.3626732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6503705) q[0];
sx q[0];
rz(-3.1351884) q[0];
sx q[0];
rz(0.48625913) q[0];
rz(2.4468415) q[1];
sx q[1];
rz(-2.7943352) q[1];
sx q[1];
rz(2.6659226) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.75695) q[0];
sx q[0];
rz(-1.6044549) q[0];
sx q[0];
rz(1.5348987) q[0];
x q[1];
rz(0.86464244) q[2];
sx q[2];
rz(-1.5448031) q[2];
sx q[2];
rz(-0.038250462) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1022288) q[1];
sx q[1];
rz(-2.1783531) q[1];
sx q[1];
rz(-1.5629014) q[1];
rz(-pi) q[2];
rz(2.6127897) q[3];
sx q[3];
rz(-0.35548726) q[3];
sx q[3];
rz(-0.58099174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4645369) q[2];
sx q[2];
rz(-3.0812283) q[2];
sx q[2];
rz(1.8439058) q[2];
rz(-1.248598) q[3];
sx q[3];
rz(-0.55515754) q[3];
sx q[3];
rz(2.7738074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1260592) q[0];
sx q[0];
rz(-1.561469) q[0];
sx q[0];
rz(-1.4373686) q[0];
rz(2.2592648) q[1];
sx q[1];
rz(-0.066451646) q[1];
sx q[1];
rz(0.8393504) q[1];
rz(2.7693314) q[2];
sx q[2];
rz(-0.098975565) q[2];
sx q[2];
rz(3.0437058) q[2];
rz(0.23918693) q[3];
sx q[3];
rz(-1.5755972) q[3];
sx q[3];
rz(1.5920873) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
