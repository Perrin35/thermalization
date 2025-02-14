OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.68706694) q[0];
sx q[0];
rz(-1.878976) q[0];
sx q[0];
rz(-0.73448056) q[0];
rz(1.6425411) q[1];
sx q[1];
rz(-1.2999111) q[1];
sx q[1];
rz(-2.1631961) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38815094) q[0];
sx q[0];
rz(-1.7072258) q[0];
sx q[0];
rz(0.0078972422) q[0];
rz(-pi) q[1];
rz(1.0389655) q[2];
sx q[2];
rz(-2.2072993) q[2];
sx q[2];
rz(-2.7920565) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.96127993) q[1];
sx q[1];
rz(-0.35113564) q[1];
sx q[1];
rz(1.9957955) q[1];
rz(1.0722854) q[3];
sx q[3];
rz(-1.6408593) q[3];
sx q[3];
rz(0.69341502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3413099) q[2];
sx q[2];
rz(-1.6494696) q[2];
sx q[2];
rz(0.71551234) q[2];
rz(-1.4095151) q[3];
sx q[3];
rz(-0.87604299) q[3];
sx q[3];
rz(1.6375665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.3187123) q[0];
sx q[0];
rz(-0.57682288) q[0];
sx q[0];
rz(3.1067644) q[0];
rz(-2.4233129) q[1];
sx q[1];
rz(-1.4052582) q[1];
sx q[1];
rz(-1.9504331) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1930192) q[0];
sx q[0];
rz(-2.5002242) q[0];
sx q[0];
rz(-1.3032342) q[0];
x q[1];
rz(-1.4227563) q[2];
sx q[2];
rz(-0.89385539) q[2];
sx q[2];
rz(-1.1803152) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9893966) q[1];
sx q[1];
rz(-2.4862444) q[1];
sx q[1];
rz(2.6776777) q[1];
x q[2];
rz(2.9510849) q[3];
sx q[3];
rz(-1.1393875) q[3];
sx q[3];
rz(2.0634758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3388153) q[2];
sx q[2];
rz(-1.85314) q[2];
sx q[2];
rz(-0.046099376) q[2];
rz(0.80373803) q[3];
sx q[3];
rz(-1.2396783) q[3];
sx q[3];
rz(0.55155915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7130724) q[0];
sx q[0];
rz(-1.5621194) q[0];
sx q[0];
rz(2.5285316) q[0];
rz(-2.8763981) q[1];
sx q[1];
rz(-1.801633) q[1];
sx q[1];
rz(-0.048096098) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4652258) q[0];
sx q[0];
rz(-2.6597659) q[0];
sx q[0];
rz(-1.5761907) q[0];
x q[1];
rz(-2.5597495) q[2];
sx q[2];
rz(-1.1238691) q[2];
sx q[2];
rz(-2.4012411) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.16471699) q[1];
sx q[1];
rz(-1.6304468) q[1];
sx q[1];
rz(3.0890134) q[1];
rz(-pi) q[2];
rz(-2.7907324) q[3];
sx q[3];
rz(-1.8985629) q[3];
sx q[3];
rz(1.3223437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4464104) q[2];
sx q[2];
rz(-1.1626652) q[2];
sx q[2];
rz(0.41333684) q[2];
rz(0.6684331) q[3];
sx q[3];
rz(-1.3276525) q[3];
sx q[3];
rz(1.7641164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17871019) q[0];
sx q[0];
rz(-0.2186192) q[0];
sx q[0];
rz(0.14337732) q[0];
rz(-0.42477056) q[1];
sx q[1];
rz(-1.2062585) q[1];
sx q[1];
rz(0.43407789) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9044943) q[0];
sx q[0];
rz(-1.6088271) q[0];
sx q[0];
rz(1.664094) q[0];
rz(2.4894478) q[2];
sx q[2];
rz(-1.4776448) q[2];
sx q[2];
rz(0.31983122) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.078316704) q[1];
sx q[1];
rz(-0.5393636) q[1];
sx q[1];
rz(-1.2961948) q[1];
rz(1.3531209) q[3];
sx q[3];
rz(-0.77131144) q[3];
sx q[3];
rz(1.2319437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.268959) q[2];
sx q[2];
rz(-2.3838398) q[2];
sx q[2];
rz(-2.7244275) q[2];
rz(-0.68552351) q[3];
sx q[3];
rz(-1.7020107) q[3];
sx q[3];
rz(-3.0912257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3454469) q[0];
sx q[0];
rz(-0.31679994) q[0];
sx q[0];
rz(1.5078804) q[0];
rz(2.4729074) q[1];
sx q[1];
rz(-1.94328) q[1];
sx q[1];
rz(1.1315469) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39172669) q[0];
sx q[0];
rz(-1.6035115) q[0];
sx q[0];
rz(-0.97264887) q[0];
rz(-2.8077233) q[2];
sx q[2];
rz(-1.0354932) q[2];
sx q[2];
rz(1.6795422) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2315438) q[1];
sx q[1];
rz(-2.0374415) q[1];
sx q[1];
rz(1.3620767) q[1];
x q[2];
rz(-0.82238166) q[3];
sx q[3];
rz(-2.9677792) q[3];
sx q[3];
rz(-1.2863152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.23182997) q[2];
sx q[2];
rz(-0.54002395) q[2];
sx q[2];
rz(2.6798999) q[2];
rz(-2.8716904) q[3];
sx q[3];
rz(-1.2609127) q[3];
sx q[3];
rz(2.4902978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5676024) q[0];
sx q[0];
rz(-0.42581588) q[0];
sx q[0];
rz(-2.0800166) q[0];
rz(0.65879446) q[1];
sx q[1];
rz(-1.7196722) q[1];
sx q[1];
rz(-0.40491358) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39570582) q[0];
sx q[0];
rz(-0.86942023) q[0];
sx q[0];
rz(1.7084136) q[0];
x q[1];
rz(1.986386) q[2];
sx q[2];
rz(-2.0743557) q[2];
sx q[2];
rz(2.1301309) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4872849) q[1];
sx q[1];
rz(-0.52738076) q[1];
sx q[1];
rz(1.2913778) q[1];
x q[2];
rz(-0.93308927) q[3];
sx q[3];
rz(-1.8130867) q[3];
sx q[3];
rz(1.0090431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0460661) q[2];
sx q[2];
rz(-2.1370856) q[2];
sx q[2];
rz(-1.3859762) q[2];
rz(2.4447377) q[3];
sx q[3];
rz(-0.98074061) q[3];
sx q[3];
rz(2.3314886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5974469) q[0];
sx q[0];
rz(-2.5283041) q[0];
sx q[0];
rz(-0.59462732) q[0];
rz(2.0095297) q[1];
sx q[1];
rz(-1.7375528) q[1];
sx q[1];
rz(-2.6304257) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88379556) q[0];
sx q[0];
rz(-0.57475677) q[0];
sx q[0];
rz(0.099959838) q[0];
rz(-pi) q[1];
rz(-2.862045) q[2];
sx q[2];
rz(-1.8862533) q[2];
sx q[2];
rz(0.73260288) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0446544) q[1];
sx q[1];
rz(-1.4739828) q[1];
sx q[1];
rz(-1.3083878) q[1];
rz(-2.0603544) q[3];
sx q[3];
rz(-0.900002) q[3];
sx q[3];
rz(1.9157992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1749997) q[2];
sx q[2];
rz(-0.24093691) q[2];
sx q[2];
rz(-0.28755507) q[2];
rz(1.1047085) q[3];
sx q[3];
rz(-1.46773) q[3];
sx q[3];
rz(-0.12565676) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6519046) q[0];
sx q[0];
rz(-0.37864417) q[0];
sx q[0];
rz(2.4878159) q[0];
rz(2.6405624) q[1];
sx q[1];
rz(-0.72622314) q[1];
sx q[1];
rz(-0.78530606) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5637999) q[0];
sx q[0];
rz(-2.1350088) q[0];
sx q[0];
rz(-1.7336506) q[0];
x q[1];
rz(-2.5748524) q[2];
sx q[2];
rz(-1.1840242) q[2];
sx q[2];
rz(-0.45505986) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.056902) q[1];
sx q[1];
rz(-0.71699079) q[1];
sx q[1];
rz(-2.8193974) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6527965) q[3];
sx q[3];
rz(-1.3424216) q[3];
sx q[3];
rz(2.5395951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.952897) q[2];
sx q[2];
rz(-1.5263824) q[2];
sx q[2];
rz(0.3696028) q[2];
rz(-2.8772307) q[3];
sx q[3];
rz(-1.2814458) q[3];
sx q[3];
rz(-0.00082357426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13067506) q[0];
sx q[0];
rz(-2.4322746) q[0];
sx q[0];
rz(1.5189019) q[0];
rz(-1.2394637) q[1];
sx q[1];
rz(-1.112097) q[1];
sx q[1];
rz(-2.1942203) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6560871) q[0];
sx q[0];
rz(-1.5099021) q[0];
sx q[0];
rz(-1.1709309) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9765509) q[2];
sx q[2];
rz(-2.7758935) q[2];
sx q[2];
rz(2.5132881) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4650146) q[1];
sx q[1];
rz(-1.7316707) q[1];
sx q[1];
rz(-1.7382998) q[1];
rz(-pi) q[2];
rz(1.4072284) q[3];
sx q[3];
rz(-1.6298098) q[3];
sx q[3];
rz(-1.7342156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6685247) q[2];
sx q[2];
rz(-0.99978414) q[2];
sx q[2];
rz(2.0295985) q[2];
rz(-0.2019349) q[3];
sx q[3];
rz(-0.58378059) q[3];
sx q[3];
rz(0.37566617) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19875459) q[0];
sx q[0];
rz(-0.18809479) q[0];
sx q[0];
rz(1.6756219) q[0];
rz(-1.5415883) q[1];
sx q[1];
rz(-1.8406638) q[1];
sx q[1];
rz(-0.25513729) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32624558) q[0];
sx q[0];
rz(-1.4370942) q[0];
sx q[0];
rz(0.017450138) q[0];
x q[1];
rz(2.3681247) q[2];
sx q[2];
rz(-1.8102526) q[2];
sx q[2];
rz(-2.4635893) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0855116) q[1];
sx q[1];
rz(-1.2930808) q[1];
sx q[1];
rz(3.059832) q[1];
rz(1.7621683) q[3];
sx q[3];
rz(-0.81052033) q[3];
sx q[3];
rz(0.22374353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.61350358) q[2];
sx q[2];
rz(-1.9518423) q[2];
sx q[2];
rz(-0.36821723) q[2];
rz(0.28299371) q[3];
sx q[3];
rz(-2.8734983) q[3];
sx q[3];
rz(2.8158269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.014864347) q[0];
sx q[0];
rz(-0.8664425) q[0];
sx q[0];
rz(2.0436825) q[0];
rz(2.6425843) q[1];
sx q[1];
rz(-1.9022763) q[1];
sx q[1];
rz(0.3872445) q[1];
rz(1.4426184) q[2];
sx q[2];
rz(-1.2405996) q[2];
sx q[2];
rz(2.2157359) q[2];
rz(-1.7718689) q[3];
sx q[3];
rz(-1.6314421) q[3];
sx q[3];
rz(-0.31753123) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
