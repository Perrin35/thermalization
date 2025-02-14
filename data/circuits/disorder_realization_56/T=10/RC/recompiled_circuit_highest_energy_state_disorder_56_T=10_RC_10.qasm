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
rz(2.7409878) q[0];
sx q[0];
rz(-0.40979835) q[0];
sx q[0];
rz(1.0709437) q[0];
rz(0.61869705) q[1];
sx q[1];
rz(-2.5735811) q[1];
sx q[1];
rz(-0.84363371) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6370359) q[0];
sx q[0];
rz(-1.8352274) q[0];
sx q[0];
rz(-0.030585551) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8461692) q[2];
sx q[2];
rz(-1.43338) q[2];
sx q[2];
rz(0.58770056) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.93701332) q[1];
sx q[1];
rz(-0.76256231) q[1];
sx q[1];
rz(-2.049905) q[1];
rz(2.0655965) q[3];
sx q[3];
rz(-1.9952979) q[3];
sx q[3];
rz(-0.44874661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7817276) q[2];
sx q[2];
rz(-2.3024776) q[2];
sx q[2];
rz(-3.1245933) q[2];
rz(-1.4012236) q[3];
sx q[3];
rz(-2.5798116) q[3];
sx q[3];
rz(1.2167654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.568999) q[0];
sx q[0];
rz(-1.8091135) q[0];
sx q[0];
rz(0.78729415) q[0];
rz(-0.17838082) q[1];
sx q[1];
rz(-1.9352103) q[1];
sx q[1];
rz(-1.9040727) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76233076) q[0];
sx q[0];
rz(-1.6325765) q[0];
sx q[0];
rz(0.5835377) q[0];
x q[1];
rz(-1.167088) q[2];
sx q[2];
rz(-2.4879527) q[2];
sx q[2];
rz(2.6622651) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.54531714) q[1];
sx q[1];
rz(-1.8870237) q[1];
sx q[1];
rz(-2.0278707) q[1];
rz(-pi) q[2];
rz(-1.848189) q[3];
sx q[3];
rz(-1.1451058) q[3];
sx q[3];
rz(-2.7077146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7201207) q[2];
sx q[2];
rz(-1.4482435) q[2];
sx q[2];
rz(-0.063701542) q[2];
rz(-1.2740159) q[3];
sx q[3];
rz(-2.2985022) q[3];
sx q[3];
rz(1.5773076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.098123) q[0];
sx q[0];
rz(-1.948057) q[0];
sx q[0];
rz(-0.7435588) q[0];
rz(2.6877563) q[1];
sx q[1];
rz(-1.791714) q[1];
sx q[1];
rz(2.7226864) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4672353) q[0];
sx q[0];
rz(-1.7667223) q[0];
sx q[0];
rz(1.118994) q[0];
rz(-pi) q[1];
rz(-1.9100283) q[2];
sx q[2];
rz(-2.3776109) q[2];
sx q[2];
rz(-0.52896777) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3943399) q[1];
sx q[1];
rz(-1.2868007) q[1];
sx q[1];
rz(-2.2153077) q[1];
x q[2];
rz(-1.4019074) q[3];
sx q[3];
rz(-0.68615507) q[3];
sx q[3];
rz(-2.9952733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.087674) q[2];
sx q[2];
rz(-2.3232465) q[2];
sx q[2];
rz(-1.9107001) q[2];
rz(2.6026717) q[3];
sx q[3];
rz(-1.6659707) q[3];
sx q[3];
rz(-1.4628791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26381275) q[0];
sx q[0];
rz(-0.75436622) q[0];
sx q[0];
rz(0.60212773) q[0];
rz(1.6944132) q[1];
sx q[1];
rz(-2.4592631) q[1];
sx q[1];
rz(0.86300659) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38463889) q[0];
sx q[0];
rz(-0.85345972) q[0];
sx q[0];
rz(1.4199281) q[0];
rz(-pi) q[1];
x q[1];
rz(2.237488) q[2];
sx q[2];
rz(-2.6764884) q[2];
sx q[2];
rz(2.1423774) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3366886) q[1];
sx q[1];
rz(-0.64283481) q[1];
sx q[1];
rz(-1.1240608) q[1];
rz(-pi) q[2];
x q[2];
rz(1.225215) q[3];
sx q[3];
rz(-0.76831078) q[3];
sx q[3];
rz(0.5925203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.065923125) q[2];
sx q[2];
rz(-1.9663591) q[2];
sx q[2];
rz(-2.2389331) q[2];
rz(0.81501189) q[3];
sx q[3];
rz(-0.60324001) q[3];
sx q[3];
rz(-0.41316113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.0937423) q[0];
sx q[0];
rz(-3.0939026) q[0];
sx q[0];
rz(-0.93511859) q[0];
rz(-1.6085666) q[1];
sx q[1];
rz(-0.32564274) q[1];
sx q[1];
rz(-0.26587048) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1810721) q[0];
sx q[0];
rz(-1.4645394) q[0];
sx q[0];
rz(1.8828859) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8498672) q[2];
sx q[2];
rz(-0.84013591) q[2];
sx q[2];
rz(2.1389291) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2468202) q[1];
sx q[1];
rz(-0.38684626) q[1];
sx q[1];
rz(0.43718613) q[1];
rz(1.3060027) q[3];
sx q[3];
rz(-1.6700498) q[3];
sx q[3];
rz(-2.3276727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.68044668) q[2];
sx q[2];
rz(-1.5004044) q[2];
sx q[2];
rz(-0.45073304) q[2];
rz(-1.9510673) q[3];
sx q[3];
rz(-0.7014941) q[3];
sx q[3];
rz(-2.7019971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9615237) q[0];
sx q[0];
rz(-0.30338767) q[0];
sx q[0];
rz(-3.0971089) q[0];
rz(-0.26875177) q[1];
sx q[1];
rz(-0.93695295) q[1];
sx q[1];
rz(-0.43089795) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6857032) q[0];
sx q[0];
rz(-2.1834032) q[0];
sx q[0];
rz(1.8926748) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0111647) q[2];
sx q[2];
rz(-0.75087386) q[2];
sx q[2];
rz(-0.73926413) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.84316501) q[1];
sx q[1];
rz(-0.80363217) q[1];
sx q[1];
rz(2.1433565) q[1];
rz(-pi) q[2];
rz(-1.5765801) q[3];
sx q[3];
rz(-0.91051379) q[3];
sx q[3];
rz(-0.74951142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0253133) q[2];
sx q[2];
rz(-0.37386027) q[2];
sx q[2];
rz(-0.6244134) q[2];
rz(2.0404909) q[3];
sx q[3];
rz(-0.81239429) q[3];
sx q[3];
rz(-0.98744923) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91200149) q[0];
sx q[0];
rz(-1.0593375) q[0];
sx q[0];
rz(0.68369317) q[0];
rz(-1.6992441) q[1];
sx q[1];
rz(-0.78386274) q[1];
sx q[1];
rz(-2.9396465) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6091222) q[0];
sx q[0];
rz(-1.2426071) q[0];
sx q[0];
rz(-3.1402324) q[0];
rz(-pi) q[1];
rz(1.5187289) q[2];
sx q[2];
rz(-2.3935742) q[2];
sx q[2];
rz(1.5128653) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8117042) q[1];
sx q[1];
rz(-1.0502019) q[1];
sx q[1];
rz(1.4938526) q[1];
rz(-pi) q[2];
rz(-0.087183909) q[3];
sx q[3];
rz(-2.6823061) q[3];
sx q[3];
rz(-2.313314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.086229) q[2];
sx q[2];
rz(-1.7217041) q[2];
sx q[2];
rz(1.8602547) q[2];
rz(2.4751439) q[3];
sx q[3];
rz(-1.0969578) q[3];
sx q[3];
rz(2.5950529) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1619038) q[0];
sx q[0];
rz(-0.55210102) q[0];
sx q[0];
rz(0.51396489) q[0];
rz(-0.95298302) q[1];
sx q[1];
rz(-0.62478462) q[1];
sx q[1];
rz(-2.0228588) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0054487451) q[0];
sx q[0];
rz(-1.4892254) q[0];
sx q[0];
rz(0.42382999) q[0];
rz(-0.4301901) q[2];
sx q[2];
rz(-0.83900827) q[2];
sx q[2];
rz(0.48626394) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8427022) q[1];
sx q[1];
rz(-1.3566847) q[1];
sx q[1];
rz(-1.9886809) q[1];
rz(-pi) q[2];
rz(-2.5932157) q[3];
sx q[3];
rz(-1.8781239) q[3];
sx q[3];
rz(0.69068324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3861367) q[2];
sx q[2];
rz(-0.35217199) q[2];
sx q[2];
rz(0.90292162) q[2];
rz(-1.4593982) q[3];
sx q[3];
rz(-1.3420339) q[3];
sx q[3];
rz(-0.82227796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4626386) q[0];
sx q[0];
rz(-2.8215388) q[0];
sx q[0];
rz(-1.1902887) q[0];
rz(-0.32084385) q[1];
sx q[1];
rz(-1.4535934) q[1];
sx q[1];
rz(-1.2215337) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4325334) q[0];
sx q[0];
rz(-1.7072868) q[0];
sx q[0];
rz(-0.20705072) q[0];
x q[1];
rz(0.83619976) q[2];
sx q[2];
rz(-1.778002) q[2];
sx q[2];
rz(1.5839603) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1561191) q[1];
sx q[1];
rz(-1.3669984) q[1];
sx q[1];
rz(2.3812095) q[1];
rz(-pi) q[2];
rz(2.9469195) q[3];
sx q[3];
rz(-1.750573) q[3];
sx q[3];
rz(0.67357066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.90423501) q[2];
sx q[2];
rz(-0.15494896) q[2];
sx q[2];
rz(-2.8083189) q[2];
rz(-2.1244369) q[3];
sx q[3];
rz(-2.1867496) q[3];
sx q[3];
rz(-2.8216383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79621133) q[0];
sx q[0];
rz(-2.6295202) q[0];
sx q[0];
rz(0.91878015) q[0];
rz(-2.9504919) q[1];
sx q[1];
rz(-2.7156576) q[1];
sx q[1];
rz(2.3474615) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0427754) q[0];
sx q[0];
rz(-2.2198027) q[0];
sx q[0];
rz(0.88943435) q[0];
rz(-2.7808519) q[2];
sx q[2];
rz(-1.7684002) q[2];
sx q[2];
rz(-0.2195356) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0550784) q[1];
sx q[1];
rz(-0.76755133) q[1];
sx q[1];
rz(1.1805736) q[1];
rz(-pi) q[2];
rz(0.32148949) q[3];
sx q[3];
rz(-1.3703385) q[3];
sx q[3];
rz(3.0908302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.11067757) q[2];
sx q[2];
rz(-2.5258749) q[2];
sx q[2];
rz(-2.3827629) q[2];
rz(-2.2714254) q[3];
sx q[3];
rz(-2.3535959) q[3];
sx q[3];
rz(-2.0738535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(3.1398685) q[0];
sx q[0];
rz(-1.4586466) q[0];
sx q[0];
rz(1.9139105) q[0];
rz(-0.43791804) q[1];
sx q[1];
rz(-1.7057849) q[1];
sx q[1];
rz(1.5954856) q[1];
rz(0.48702892) q[2];
sx q[2];
rz(-1.9430046) q[2];
sx q[2];
rz(-0.76443048) q[2];
rz(-2.6043456) q[3];
sx q[3];
rz(-0.8575079) q[3];
sx q[3];
rz(-1.433123) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
