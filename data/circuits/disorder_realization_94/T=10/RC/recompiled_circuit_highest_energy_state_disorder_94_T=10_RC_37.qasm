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
rz(1.2430159) q[0];
sx q[0];
rz(-1.1216811) q[0];
sx q[0];
rz(-2.6774874) q[0];
rz(2.3857181) q[1];
sx q[1];
rz(2.2321489) q[1];
sx q[1];
rz(7.5997054) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6073734) q[0];
sx q[0];
rz(-1.6053622) q[0];
sx q[0];
rz(3.0764293) q[0];
rz(-pi) q[1];
rz(1.520697) q[2];
sx q[2];
rz(-1.9110137) q[2];
sx q[2];
rz(-2.9638247) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2006372) q[1];
sx q[1];
rz(-0.66087729) q[1];
sx q[1];
rz(0.4502859) q[1];
x q[2];
rz(-1.2795696) q[3];
sx q[3];
rz(-1.816245) q[3];
sx q[3];
rz(-0.042118532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.89062771) q[2];
sx q[2];
rz(-0.35718063) q[2];
sx q[2];
rz(2.8446021) q[2];
rz(-0.65303981) q[3];
sx q[3];
rz(-1.5831455) q[3];
sx q[3];
rz(-2.0853341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62946573) q[0];
sx q[0];
rz(-1.0221721) q[0];
sx q[0];
rz(-0.28489354) q[0];
rz(-3.0902872) q[1];
sx q[1];
rz(-2.3785794) q[1];
sx q[1];
rz(0.6368534) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013971) q[0];
sx q[0];
rz(-0.68249615) q[0];
sx q[0];
rz(-2.3861814) q[0];
rz(-1.571973) q[2];
sx q[2];
rz(-1.92244) q[2];
sx q[2];
rz(-0.52368173) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0242549) q[1];
sx q[1];
rz(-1.7452612) q[1];
sx q[1];
rz(-1.6478478) q[1];
x q[2];
rz(-0.15588197) q[3];
sx q[3];
rz(-0.71195275) q[3];
sx q[3];
rz(-0.32526325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.18179831) q[2];
sx q[2];
rz(-0.85796285) q[2];
sx q[2];
rz(-2.174343) q[2];
rz(-0.26425427) q[3];
sx q[3];
rz(-1.6382917) q[3];
sx q[3];
rz(1.2359469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13977519) q[0];
sx q[0];
rz(-1.5793261) q[0];
sx q[0];
rz(-1.1199957) q[0];
rz(-2.6657875) q[1];
sx q[1];
rz(-2.0640524) q[1];
sx q[1];
rz(-0.30398223) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8037) q[0];
sx q[0];
rz(-2.5753382) q[0];
sx q[0];
rz(-1.5192274) q[0];
rz(-1.4690222) q[2];
sx q[2];
rz(-1.3341122) q[2];
sx q[2];
rz(-2.9477811) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.099977255) q[1];
sx q[1];
rz(-1.4681889) q[1];
sx q[1];
rz(0.85081886) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.68328698) q[3];
sx q[3];
rz(-1.4158472) q[3];
sx q[3];
rz(0.72282997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1659282) q[2];
sx q[2];
rz(-2.1037481) q[2];
sx q[2];
rz(-2.4456639) q[2];
rz(1.1337229) q[3];
sx q[3];
rz(-1.1997831) q[3];
sx q[3];
rz(0.55293647) q[3];
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
rz(pi/2) q[3];
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
rz(0.43941471) q[0];
sx q[0];
rz(-2.1372097) q[0];
sx q[0];
rz(1.0473921) q[0];
rz(-1.9678496) q[1];
sx q[1];
rz(-0.67432299) q[1];
sx q[1];
rz(-1.5210927) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5562963) q[0];
sx q[0];
rz(-1.4504787) q[0];
sx q[0];
rz(-0.5043982) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9416472) q[2];
sx q[2];
rz(-2.4003138) q[2];
sx q[2];
rz(1.0560869) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.92482216) q[1];
sx q[1];
rz(-2.3433279) q[1];
sx q[1];
rz(-0.63191788) q[1];
rz(-pi) q[2];
rz(-1.2866045) q[3];
sx q[3];
rz(-1.2582558) q[3];
sx q[3];
rz(-0.72060637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0461222) q[2];
sx q[2];
rz(-1.0669402) q[2];
sx q[2];
rz(-0.53786892) q[2];
rz(1.6543903) q[3];
sx q[3];
rz(-0.40707773) q[3];
sx q[3];
rz(2.4552104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2122413) q[0];
sx q[0];
rz(-0.23433267) q[0];
sx q[0];
rz(0.60337639) q[0];
rz(0.76250184) q[1];
sx q[1];
rz(-1.1926032) q[1];
sx q[1];
rz(-2.4286483) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.514587) q[0];
sx q[0];
rz(-2.7381185) q[0];
sx q[0];
rz(1.7395942) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1080526) q[2];
sx q[2];
rz(-1.4456835) q[2];
sx q[2];
rz(2.3274328) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0281369) q[1];
sx q[1];
rz(-2.0638247) q[1];
sx q[1];
rz(-0.22844577) q[1];
rz(0.98865786) q[3];
sx q[3];
rz(-2.1063559) q[3];
sx q[3];
rz(2.5861507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.9731628) q[2];
sx q[2];
rz(-0.91732401) q[2];
sx q[2];
rz(1.6707576) q[2];
rz(2.6245608) q[3];
sx q[3];
rz(-1.0443338) q[3];
sx q[3];
rz(1.8487336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9789199) q[0];
sx q[0];
rz(-2.6277268) q[0];
sx q[0];
rz(-0.55147076) q[0];
rz(-3.1061213) q[1];
sx q[1];
rz(-2.008581) q[1];
sx q[1];
rz(1.5303401) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96824232) q[0];
sx q[0];
rz(-2.1538072) q[0];
sx q[0];
rz(-0.75847404) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6257203) q[2];
sx q[2];
rz(-0.74797219) q[2];
sx q[2];
rz(-2.2750281) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2179709) q[1];
sx q[1];
rz(-1.21232) q[1];
sx q[1];
rz(2.3450407) q[1];
rz(-pi) q[2];
rz(-0.91776589) q[3];
sx q[3];
rz(-1.5035331) q[3];
sx q[3];
rz(1.3239087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0643206) q[2];
sx q[2];
rz(-2.6214226) q[2];
sx q[2];
rz(-1.2426097) q[2];
rz(0.51314917) q[3];
sx q[3];
rz(-2.7276701) q[3];
sx q[3];
rz(-1.7665524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4921017) q[0];
sx q[0];
rz(-2.212337) q[0];
sx q[0];
rz(3.0249366) q[0];
rz(0.77313441) q[1];
sx q[1];
rz(-1.9832858) q[1];
sx q[1];
rz(1.6366417) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3706544) q[0];
sx q[0];
rz(-2.2794834) q[0];
sx q[0];
rz(0.38972008) q[0];
x q[1];
rz(-1.9595615) q[2];
sx q[2];
rz(-1.5156399) q[2];
sx q[2];
rz(0.15832947) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4915858) q[1];
sx q[1];
rz(-1.9284231) q[1];
sx q[1];
rz(-2.3050344) q[1];
rz(-pi) q[2];
rz(-1.2833474) q[3];
sx q[3];
rz(-2.2032167) q[3];
sx q[3];
rz(0.27770628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2152805) q[2];
sx q[2];
rz(-0.24682385) q[2];
sx q[2];
rz(2.2543294) q[2];
rz(-0.67982802) q[3];
sx q[3];
rz(-2.4726548) q[3];
sx q[3];
rz(0.63290709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.585084) q[0];
sx q[0];
rz(-1.8483138) q[0];
sx q[0];
rz(2.4853117) q[0];
rz(-0.20690021) q[1];
sx q[1];
rz(-1.921939) q[1];
sx q[1];
rz(1.2178749) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.525187) q[0];
sx q[0];
rz(-1.5357247) q[0];
sx q[0];
rz(3.0588018) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98274173) q[2];
sx q[2];
rz(-1.8652735) q[2];
sx q[2];
rz(-0.5628995) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6177442) q[1];
sx q[1];
rz(-0.19568014) q[1];
sx q[1];
rz(-0.62472549) q[1];
rz(-2.5303305) q[3];
sx q[3];
rz(-1.3484133) q[3];
sx q[3];
rz(-2.5739447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8021585) q[2];
sx q[2];
rz(-1.3863486) q[2];
sx q[2];
rz(0.68186861) q[2];
rz(-1.1784461) q[3];
sx q[3];
rz(-0.75740564) q[3];
sx q[3];
rz(-1.9044378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7449529) q[0];
sx q[0];
rz(-1.5322026) q[0];
sx q[0];
rz(-2.9314281) q[0];
rz(2.919803) q[1];
sx q[1];
rz(-2.2237325) q[1];
sx q[1];
rz(0.04235696) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6307459) q[0];
sx q[0];
rz(-1.5928942) q[0];
sx q[0];
rz(-0.18760292) q[0];
x q[1];
rz(0.63318166) q[2];
sx q[2];
rz(-2.0821619) q[2];
sx q[2];
rz(2.2610469) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6575129) q[1];
sx q[1];
rz(-1.8575605) q[1];
sx q[1];
rz(-2.5628255) q[1];
x q[2];
rz(0.70019763) q[3];
sx q[3];
rz(-2.4841016) q[3];
sx q[3];
rz(1.0669062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.28045851) q[2];
sx q[2];
rz(-1.0886085) q[2];
sx q[2];
rz(-2.4750278) q[2];
rz(2.376453) q[3];
sx q[3];
rz(-2.9355526) q[3];
sx q[3];
rz(1.4008745) q[3];
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
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53369451) q[0];
sx q[0];
rz(-0.38132897) q[0];
sx q[0];
rz(-2.4545942) q[0];
rz(-1.8580565) q[1];
sx q[1];
rz(-1.9891519) q[1];
sx q[1];
rz(2.5680465) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.682293) q[0];
sx q[0];
rz(-2.1577765) q[0];
sx q[0];
rz(-2.1989376) q[0];
rz(-1.1186662) q[2];
sx q[2];
rz(-1.3616478) q[2];
sx q[2];
rz(0.5224519) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.95038762) q[1];
sx q[1];
rz(-2.5053686) q[1];
sx q[1];
rz(2.1649182) q[1];
rz(-pi) q[2];
rz(2.2881094) q[3];
sx q[3];
rz(-0.9694582) q[3];
sx q[3];
rz(-0.36313148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0987229) q[2];
sx q[2];
rz(-1.9645773) q[2];
sx q[2];
rz(-3.0549808) q[2];
rz(-3.0232271) q[3];
sx q[3];
rz(-1.5990853) q[3];
sx q[3];
rz(-3.0983483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31556986) q[0];
sx q[0];
rz(-1.0181027) q[0];
sx q[0];
rz(-2.3636567) q[0];
rz(2.8553873) q[1];
sx q[1];
rz(-2.310391) q[1];
sx q[1];
rz(-1.8117767) q[1];
rz(0.33230843) q[2];
sx q[2];
rz(-1.341991) q[2];
sx q[2];
rz(1.1534635) q[2];
rz(0.44888857) q[3];
sx q[3];
rz(-1.8406624) q[3];
sx q[3];
rz(0.016907666) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
