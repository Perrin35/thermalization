OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0491068) q[0];
sx q[0];
rz(4.4409039) q[0];
sx q[0];
rz(9.470603) q[0];
rz(-1.8419645) q[1];
sx q[1];
rz(-1.3388495) q[1];
sx q[1];
rz(-1.2531228) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078431167) q[0];
sx q[0];
rz(-0.51126152) q[0];
sx q[0];
rz(-1.0529165) q[0];
x q[1];
rz(-1.5270751) q[2];
sx q[2];
rz(-1.3695903) q[2];
sx q[2];
rz(-1.218856) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.33838446) q[1];
sx q[1];
rz(-0.82370355) q[1];
sx q[1];
rz(1.7335152) q[1];
x q[2];
rz(-1.6982444) q[3];
sx q[3];
rz(-1.5747254) q[3];
sx q[3];
rz(-1.9779825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.32156285) q[2];
sx q[2];
rz(-1.7367312) q[2];
sx q[2];
rz(-2.826214) q[2];
rz(-0.086221181) q[3];
sx q[3];
rz(-3.0486221) q[3];
sx q[3];
rz(2.2975547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4013937) q[0];
sx q[0];
rz(-1.711015) q[0];
sx q[0];
rz(0.33541086) q[0];
rz(-0.3849349) q[1];
sx q[1];
rz(-0.79610577) q[1];
sx q[1];
rz(-1.8523432) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9612564) q[0];
sx q[0];
rz(-1.633005) q[0];
sx q[0];
rz(0.38952413) q[0];
rz(1.4440437) q[2];
sx q[2];
rz(-2.4962262) q[2];
sx q[2];
rz(2.5996948) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.34270479) q[1];
sx q[1];
rz(-1.2884166) q[1];
sx q[1];
rz(-2.7434231) q[1];
x q[2];
rz(-0.13522526) q[3];
sx q[3];
rz(-1.8117419) q[3];
sx q[3];
rz(-1.8640445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5424767) q[2];
sx q[2];
rz(-1.3834388) q[2];
sx q[2];
rz(-1.6790338) q[2];
rz(2.2595432) q[3];
sx q[3];
rz(-0.65989152) q[3];
sx q[3];
rz(-1.6409469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7733234) q[0];
sx q[0];
rz(-0.48650807) q[0];
sx q[0];
rz(2.4988556) q[0];
rz(2.2679988) q[1];
sx q[1];
rz(-1.830955) q[1];
sx q[1];
rz(2.1547623) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1571647) q[0];
sx q[0];
rz(-1.720795) q[0];
sx q[0];
rz(-2.2638129) q[0];
rz(-pi) q[1];
rz(-1.6880694) q[2];
sx q[2];
rz(-2.6542695) q[2];
sx q[2];
rz(2.8685958) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.22443188) q[1];
sx q[1];
rz(-0.81714918) q[1];
sx q[1];
rz(2.3520873) q[1];
rz(1.945472) q[3];
sx q[3];
rz(-1.5244532) q[3];
sx q[3];
rz(-2.1425193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.70283908) q[2];
sx q[2];
rz(-3.1162016) q[2];
sx q[2];
rz(1.0566443) q[2];
rz(1.3422525) q[3];
sx q[3];
rz(-1.4547576) q[3];
sx q[3];
rz(-2.2630283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.037828) q[0];
sx q[0];
rz(-0.28452888) q[0];
sx q[0];
rz(-1.9985265) q[0];
rz(0.68483886) q[1];
sx q[1];
rz(-1.3416483) q[1];
sx q[1];
rz(-0.24982223) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1764602) q[0];
sx q[0];
rz(-1.1145381) q[0];
sx q[0];
rz(0.81911032) q[0];
rz(-pi) q[1];
rz(-2.5673713) q[2];
sx q[2];
rz(-1.5674233) q[2];
sx q[2];
rz(-0.33004883) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.22801655) q[1];
sx q[1];
rz(-0.98015112) q[1];
sx q[1];
rz(0.98286079) q[1];
rz(-3.1400348) q[3];
sx q[3];
rz(-1.8127725) q[3];
sx q[3];
rz(2.2037857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7217241) q[2];
sx q[2];
rz(-1.7896174) q[2];
sx q[2];
rz(-1.295759) q[2];
rz(1.7660247) q[3];
sx q[3];
rz(-1.7121168) q[3];
sx q[3];
rz(3.0425369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7527723) q[0];
sx q[0];
rz(-1.3714014) q[0];
sx q[0];
rz(-0.8465299) q[0];
rz(1.9580152) q[1];
sx q[1];
rz(-1.1776244) q[1];
sx q[1];
rz(1.902098) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1355818) q[0];
sx q[0];
rz(-2.860489) q[0];
sx q[0];
rz(1.1233888) q[0];
rz(-pi) q[1];
rz(-1.3888277) q[2];
sx q[2];
rz(-1.2335868) q[2];
sx q[2];
rz(0.47949258) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2927276) q[1];
sx q[1];
rz(-1.5827109) q[1];
sx q[1];
rz(-1.6267651) q[1];
rz(-pi) q[2];
rz(-1.4828478) q[3];
sx q[3];
rz(-2.5426358) q[3];
sx q[3];
rz(0.14725895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6401297) q[2];
sx q[2];
rz(-1.632246) q[2];
sx q[2];
rz(-2.787369) q[2];
rz(2.4223478) q[3];
sx q[3];
rz(-2.5651599) q[3];
sx q[3];
rz(2.332212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.7648014) q[0];
sx q[0];
rz(-2.9053595) q[0];
sx q[0];
rz(-2.7948622) q[0];
rz(2.4193343) q[1];
sx q[1];
rz(-1.6112695) q[1];
sx q[1];
rz(-2.9647656) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5520598) q[0];
sx q[0];
rz(-2.0314275) q[0];
sx q[0];
rz(-0.15710196) q[0];
rz(-pi) q[1];
rz(3.0484285) q[2];
sx q[2];
rz(-1.0736862) q[2];
sx q[2];
rz(0.99592956) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.8521148) q[1];
sx q[1];
rz(-1.7758217) q[1];
sx q[1];
rz(-1.4641028) q[1];
rz(-pi) q[2];
rz(-1.7309639) q[3];
sx q[3];
rz(-2.5729542) q[3];
sx q[3];
rz(1.7344432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5727545) q[2];
sx q[2];
rz(-0.73358959) q[2];
sx q[2];
rz(-1.0432165) q[2];
rz(0.15626945) q[3];
sx q[3];
rz(-1.0633435) q[3];
sx q[3];
rz(-0.094303057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8909376) q[0];
sx q[0];
rz(-2.1112878) q[0];
sx q[0];
rz(-1.913273) q[0];
rz(0.43371513) q[1];
sx q[1];
rz(-0.80077306) q[1];
sx q[1];
rz(2.0965651) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6390594) q[0];
sx q[0];
rz(-0.21083388) q[0];
sx q[0];
rz(1.8502214) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3190193) q[2];
sx q[2];
rz(-1.9299091) q[2];
sx q[2];
rz(2.1445779) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4853128) q[1];
sx q[1];
rz(-1.4569786) q[1];
sx q[1];
rz(0.95358221) q[1];
rz(0.52900903) q[3];
sx q[3];
rz(-2.2914791) q[3];
sx q[3];
rz(1.0320272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.23994437) q[2];
sx q[2];
rz(-2.1176691) q[2];
sx q[2];
rz(0.18590064) q[2];
rz(1.9013885) q[3];
sx q[3];
rz(-1.600772) q[3];
sx q[3];
rz(1.0143636) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1321201) q[0];
sx q[0];
rz(-1.2484231) q[0];
sx q[0];
rz(2.9820719) q[0];
rz(-0.80687833) q[1];
sx q[1];
rz(-1.9069549) q[1];
sx q[1];
rz(3.1330915) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2293866) q[0];
sx q[0];
rz(-1.1649781) q[0];
sx q[0];
rz(-0.1677558) q[0];
rz(-pi) q[1];
x q[1];
rz(1.464017) q[2];
sx q[2];
rz(-2.3651383) q[2];
sx q[2];
rz(1.8242893) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2841683) q[1];
sx q[1];
rz(-1.4393974) q[1];
sx q[1];
rz(-0.22716503) q[1];
x q[2];
rz(1.0866585) q[3];
sx q[3];
rz(-1.414742) q[3];
sx q[3];
rz(2.597996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8622417) q[2];
sx q[2];
rz(-1.778435) q[2];
sx q[2];
rz(-2.9478493) q[2];
rz(-1.6019609) q[3];
sx q[3];
rz(-1.1774457) q[3];
sx q[3];
rz(-2.0120473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2368471) q[0];
sx q[0];
rz(-1.3847677) q[0];
sx q[0];
rz(0.73614502) q[0];
rz(0.88788095) q[1];
sx q[1];
rz(-0.45279756) q[1];
sx q[1];
rz(3.091541) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5745478) q[0];
sx q[0];
rz(-2.7549681) q[0];
sx q[0];
rz(0.48166122) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29590891) q[2];
sx q[2];
rz(-0.32797932) q[2];
sx q[2];
rz(-0.28959549) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1290196) q[1];
sx q[1];
rz(-2.7677849) q[1];
sx q[1];
rz(-2.5555948) q[1];
x q[2];
rz(0.95376905) q[3];
sx q[3];
rz(-2.0514384) q[3];
sx q[3];
rz(-0.98508376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8686409) q[2];
sx q[2];
rz(-1.4401888) q[2];
sx q[2];
rz(3.1325373) q[2];
rz(-2.791259) q[3];
sx q[3];
rz(-2.83559) q[3];
sx q[3];
rz(-0.30456021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7250799) q[0];
sx q[0];
rz(-2.0820936) q[0];
sx q[0];
rz(-2.1685725) q[0];
rz(-1.5618207) q[1];
sx q[1];
rz(-0.90619722) q[1];
sx q[1];
rz(-0.80950338) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12347977) q[0];
sx q[0];
rz(-0.84313376) q[0];
sx q[0];
rz(0.48771702) q[0];
rz(-pi) q[1];
rz(0.87363957) q[2];
sx q[2];
rz(-0.29549797) q[2];
sx q[2];
rz(0.9852315) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0187518) q[1];
sx q[1];
rz(-2.0301308) q[1];
sx q[1];
rz(0.18045119) q[1];
rz(-pi) q[2];
rz(-0.9904434) q[3];
sx q[3];
rz(-1.9706315) q[3];
sx q[3];
rz(0.70580182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.96106625) q[2];
sx q[2];
rz(-0.82225353) q[2];
sx q[2];
rz(2.6194438) q[2];
rz(-0.3565878) q[3];
sx q[3];
rz(-2.5643189) q[3];
sx q[3];
rz(0.063551158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6575573) q[0];
sx q[0];
rz(-2.7900896) q[0];
sx q[0];
rz(-1.8650613) q[0];
rz(2.7784078) q[1];
sx q[1];
rz(-2.6138432) q[1];
sx q[1];
rz(1.6539727) q[1];
rz(0.54375081) q[2];
sx q[2];
rz(-2.2866572) q[2];
sx q[2];
rz(0.86845749) q[2];
rz(-1.894886) q[3];
sx q[3];
rz(-2.1264599) q[3];
sx q[3];
rz(2.651941) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
