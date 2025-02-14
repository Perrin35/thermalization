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
rz(-2.812204) q[0];
sx q[0];
rz(-0.63151276) q[0];
sx q[0];
rz(3.1407177) q[0];
rz(-0.64088351) q[1];
sx q[1];
rz(5.2823245) q[1];
sx q[1];
rz(9.7699788) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3314197) q[0];
sx q[0];
rz(-1.4799444) q[0];
sx q[0];
rz(-0.024837107) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37981793) q[2];
sx q[2];
rz(-1.4336515) q[2];
sx q[2];
rz(-1.6452546) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6908603) q[1];
sx q[1];
rz(-1.9885411) q[1];
sx q[1];
rz(-0.91207204) q[1];
rz(-pi) q[2];
rz(1.6451938) q[3];
sx q[3];
rz(-1.5545903) q[3];
sx q[3];
rz(-2.5287573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.22902809) q[2];
sx q[2];
rz(-1.3338858) q[2];
sx q[2];
rz(2.933617) q[2];
rz(0.29911706) q[3];
sx q[3];
rz(-2.5695473) q[3];
sx q[3];
rz(1.1408898) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3308554) q[0];
sx q[0];
rz(-2.9006697) q[0];
sx q[0];
rz(0.99288565) q[0];
rz(-1.8244686) q[1];
sx q[1];
rz(-2.8254852) q[1];
sx q[1];
rz(0.77450007) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82371432) q[0];
sx q[0];
rz(-1.5817989) q[0];
sx q[0];
rz(-1.7490134) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51807816) q[2];
sx q[2];
rz(-0.87730125) q[2];
sx q[2];
rz(-2.3665646) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6933813) q[1];
sx q[1];
rz(-0.44084545) q[1];
sx q[1];
rz(-2.5930659) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9162972) q[3];
sx q[3];
rz(-1.1114128) q[3];
sx q[3];
rz(0.53250203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.896686) q[2];
sx q[2];
rz(-1.2865571) q[2];
sx q[2];
rz(-2.1265105) q[2];
rz(2.8924938) q[3];
sx q[3];
rz(-2.2738012) q[3];
sx q[3];
rz(-2.7740313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86431137) q[0];
sx q[0];
rz(-2.4503777) q[0];
sx q[0];
rz(2.9840898) q[0];
rz(1.0109673) q[1];
sx q[1];
rz(-0.33282655) q[1];
sx q[1];
rz(-0.13883042) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94731733) q[0];
sx q[0];
rz(-1.9324158) q[0];
sx q[0];
rz(0.34644923) q[0];
rz(-pi) q[1];
rz(-0.46926083) q[2];
sx q[2];
rz(-2.3692694) q[2];
sx q[2];
rz(-2.0283716) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.638354) q[1];
sx q[1];
rz(-1.279083) q[1];
sx q[1];
rz(-2.9364763) q[1];
rz(-pi) q[2];
x q[2];
rz(0.54581996) q[3];
sx q[3];
rz(-0.32836093) q[3];
sx q[3];
rz(-0.92121802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.46131721) q[2];
sx q[2];
rz(-0.91247827) q[2];
sx q[2];
rz(-0.3581363) q[2];
rz(-2.4628468) q[3];
sx q[3];
rz(-0.94921422) q[3];
sx q[3];
rz(-2.1024735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(3.0963652) q[0];
sx q[0];
rz(-0.97310936) q[0];
sx q[0];
rz(2.8367693) q[0];
rz(0.40799704) q[1];
sx q[1];
rz(-1.7034737) q[1];
sx q[1];
rz(2.2136484) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.613425) q[0];
sx q[0];
rz(-1.4425462) q[0];
sx q[0];
rz(-1.3750465) q[0];
rz(-pi) q[1];
rz(-0.60073845) q[2];
sx q[2];
rz(-1.033412) q[2];
sx q[2];
rz(-1.635765) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2603307) q[1];
sx q[1];
rz(-1.7219543) q[1];
sx q[1];
rz(0.69044729) q[1];
x q[2];
rz(2.5119588) q[3];
sx q[3];
rz(-2.0009407) q[3];
sx q[3];
rz(0.54508524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1912332) q[2];
sx q[2];
rz(-2.5706036) q[2];
sx q[2];
rz(-0.18816571) q[2];
rz(1.865271) q[3];
sx q[3];
rz(-1.8264344) q[3];
sx q[3];
rz(-1.7108542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4417878) q[0];
sx q[0];
rz(-2.7975174) q[0];
sx q[0];
rz(-3.1222043) q[0];
rz(2.0806606) q[1];
sx q[1];
rz(-1.7273936) q[1];
sx q[1];
rz(-1.9020938) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0584377) q[0];
sx q[0];
rz(-0.64955901) q[0];
sx q[0];
rz(-2.044528) q[0];
x q[1];
rz(-0.6940191) q[2];
sx q[2];
rz(-2.2713714) q[2];
sx q[2];
rz(-1.4365591) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0026999) q[1];
sx q[1];
rz(-2.5475916) q[1];
sx q[1];
rz(2.6245489) q[1];
rz(-1.0723713) q[3];
sx q[3];
rz(-1.2066505) q[3];
sx q[3];
rz(-0.92723234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0197319) q[2];
sx q[2];
rz(-1.8283565) q[2];
sx q[2];
rz(-1.8149553) q[2];
rz(2.895368) q[3];
sx q[3];
rz(-1.1028057) q[3];
sx q[3];
rz(-0.8518014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5354079) q[0];
sx q[0];
rz(-2.1562205) q[0];
sx q[0];
rz(0.73915172) q[0];
rz(-0.34573653) q[1];
sx q[1];
rz(-1.3290936) q[1];
sx q[1];
rz(-2.2703222) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5369878) q[0];
sx q[0];
rz(-0.78442998) q[0];
sx q[0];
rz(-3.0434199) q[0];
rz(-0.61001444) q[2];
sx q[2];
rz(-2.2918309) q[2];
sx q[2];
rz(0.80243669) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.516219) q[1];
sx q[1];
rz(-1.1828848) q[1];
sx q[1];
rz(-0.078887786) q[1];
rz(0.63251791) q[3];
sx q[3];
rz(-2.4261203) q[3];
sx q[3];
rz(0.6882489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8295916) q[2];
sx q[2];
rz(-0.63903725) q[2];
sx q[2];
rz(2.269022) q[2];
rz(-1.568659) q[3];
sx q[3];
rz(-2.5376153) q[3];
sx q[3];
rz(2.2085371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0503814) q[0];
sx q[0];
rz(-2.8822883) q[0];
sx q[0];
rz(2.3799489) q[0];
rz(-2.7161982) q[1];
sx q[1];
rz(-1.1401221) q[1];
sx q[1];
rz(0.92686191) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2970718) q[0];
sx q[0];
rz(-1.8064587) q[0];
sx q[0];
rz(1.8051487) q[0];
x q[1];
rz(0.71304597) q[2];
sx q[2];
rz(-2.2879061) q[2];
sx q[2];
rz(0.46592679) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7898579) q[1];
sx q[1];
rz(-0.42966336) q[1];
sx q[1];
rz(-2.1273145) q[1];
rz(-pi) q[2];
rz(3.0934392) q[3];
sx q[3];
rz(-0.94814903) q[3];
sx q[3];
rz(0.66648713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0285792) q[2];
sx q[2];
rz(-2.4455652) q[2];
sx q[2];
rz(2.2704303) q[2];
rz(2.8431559) q[3];
sx q[3];
rz(-1.2454183) q[3];
sx q[3];
rz(-1.3309853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7262064) q[0];
sx q[0];
rz(-0.50006777) q[0];
sx q[0];
rz(0.4183847) q[0];
rz(2.8437974) q[1];
sx q[1];
rz(-1.9458385) q[1];
sx q[1];
rz(0.13430886) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9878085) q[0];
sx q[0];
rz(-1.0414818) q[0];
sx q[0];
rz(2.3321926) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9761132) q[2];
sx q[2];
rz(-1.3994872) q[2];
sx q[2];
rz(1.5660945) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1301535) q[1];
sx q[1];
rz(-1.7053889) q[1];
sx q[1];
rz(2.9573739) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8577096) q[3];
sx q[3];
rz(-2.0238658) q[3];
sx q[3];
rz(0.69597746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.40992752) q[2];
sx q[2];
rz(-1.8525367) q[2];
sx q[2];
rz(-0.64201391) q[2];
rz(1.0632473) q[3];
sx q[3];
rz(-2.4898873) q[3];
sx q[3];
rz(2.1003907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5409656) q[0];
sx q[0];
rz(-3.0525115) q[0];
sx q[0];
rz(-0.11216057) q[0];
rz(-1.3345831) q[1];
sx q[1];
rz(-2.5490675) q[1];
sx q[1];
rz(2.1122011) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.063113041) q[0];
sx q[0];
rz(-1.6029583) q[0];
sx q[0];
rz(-0.86634791) q[0];
rz(-1.6104944) q[2];
sx q[2];
rz(-1.2671114) q[2];
sx q[2];
rz(-3.1162709) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.16949108) q[1];
sx q[1];
rz(-1.6147366) q[1];
sx q[1];
rz(1.8010587) q[1];
rz(-pi) q[2];
rz(2.6736027) q[3];
sx q[3];
rz(-1.4889928) q[3];
sx q[3];
rz(1.8379267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.42334291) q[2];
sx q[2];
rz(-0.074904718) q[2];
sx q[2];
rz(3.1035799) q[2];
rz(-2.5293317) q[3];
sx q[3];
rz(-2.2050048) q[3];
sx q[3];
rz(3.0003701) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24899471) q[0];
sx q[0];
rz(-1.6706415) q[0];
sx q[0];
rz(0.41879642) q[0];
rz(-1.2576125) q[1];
sx q[1];
rz(-1.6323171) q[1];
sx q[1];
rz(-2.9990101) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2217956) q[0];
sx q[0];
rz(-0.16079535) q[0];
sx q[0];
rz(2.9662786) q[0];
rz(-1.293574) q[2];
sx q[2];
rz(-2.1810594) q[2];
sx q[2];
rz(2.0137613) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8293162) q[1];
sx q[1];
rz(-1.3061151) q[1];
sx q[1];
rz(2.1618202) q[1];
x q[2];
rz(2.7975818) q[3];
sx q[3];
rz(-2.0075402) q[3];
sx q[3];
rz(-0.94319447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7964145) q[2];
sx q[2];
rz(-3.0609481) q[2];
sx q[2];
rz(-0.26819116) q[2];
rz(-1.4043407) q[3];
sx q[3];
rz(-1.0920478) q[3];
sx q[3];
rz(-1.0442737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(2.0161229) q[0];
sx q[0];
rz(-2.7653427) q[0];
sx q[0];
rz(0.34633037) q[0];
rz(3.012433) q[1];
sx q[1];
rz(-1.8864514) q[1];
sx q[1];
rz(1.4056978) q[1];
rz(0.62025537) q[2];
sx q[2];
rz(-0.56369416) q[2];
sx q[2];
rz(-0.77947215) q[2];
rz(-1.2760194) q[3];
sx q[3];
rz(-1.634292) q[3];
sx q[3];
rz(2.5802142) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
